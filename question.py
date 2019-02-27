import pandas as pd
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import numpy as np
from alignment import read_in_cost_matrix,get_score_dfs

def question_one_pt_one(pd_collect):
    '''
    get which BLOSUM50 gap/ext combo run has the best false positive rate
    '''
    # get all the BLOSUM50 matrices
    blosum_fifty = pd_collect.loc[pd_collect['matrix'] == 'BLOSUM50']

    # get the index of the false positive with the best score
    max_index = blosum_fifty['false'].idxmax()

    # get the gap and extension values for that best false positive score
    # 2 for gap column, 3 for ext column
    best_gap = blosum_fifty['gap'][max_index]
    best_ext = blosum_fifty['ext'][max_index]
    return best_gap,best_ext

def plot_roc_curve(pd_df,scores,labels,name,file):
    '''
    info about understangin roc http://gim.unmc.edu/dxtests/roc3.htm
    plot receiver operating characteristic (ROC) curve for each matrix
    alternatively try https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
    '''
    sns.set_style('ticks')
    sns.set_palette("husl")

    # get the matrices
    matrix_names = pd_df['matrix'].tolist()
    unique_names = set(matrix_names)

    # for each matrix plot roc
    for m in unique_names:
        print(m)
        df = pd_df.loc[pd_df['matrix'] == m]
        fpr,tpr,threshold = roc_curve(df[labels],df[scores])
        area_under_roc_curve = auc(fpr,tpr)
        sns.lineplot(fpr,tpr,label='{0}, AUC={1}'.format(m,round(area_under_roc_curve,3)))
    plt.title('ROC curves for {0}'.format(name))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend()

    # save plot to image directory
    outdir = pathlib.Path('images')
    outfile = outdir / 'ROC_{0}.png'.format(file)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outfile),format='png')
    plt.close()

def question_one_pt_two(pd_collect,pd_scores,best_gap,best_ext):
    '''
    gap/ext values from question one pt one, which matrix has best false positive score
    '''
    # get a subset for those with the gap/ext values from the best false positive values from BLOSUM50
    best_gap_and_ext = pd_collect.loc[(pd_collect['gap']==best_gap) & (pd_collect['ext']==best_ext)]
    # get the matrix that has the best false positive rates
    max_index = best_gap_and_ext['false'].idxmax()

    # get the matrix that from that index, 4 for matrix column
    best_matrix = best_gap_and_ext['matrix'][max_index]

    plot_roc_curve(pd_scores,'scores','labels','best performing matrix, {0}'.format(best_matrix),'Best_'.format(best_matrix))

    return best_matrix

def question_one_pt_three(pd_scores,best_matrix):
    '''
    normalize smith and waterman score by length of shortest pair sequence
    roc curves for best matrix, normalized scores
    '''

    plot_roc_curve(pd_scores,'norm_scores','labels','normalized scores, {0}'.format(best_matrix),'Norm_'.format(best_matrix))

def plot_distributions(pd_df,matrix,threshold):
    '''
    plot the false and true positives around the threshold
    '''
    '''
    info about understangin roc http://gim.unmc.edu/dxtests/roc3.htm
    plot receiver operating characteristic (ROC) curve for each matrix
    alternatively try https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
    '''
    sns.set_style('ticks')
    sns.set_palette("husl")

    df = pd_df.loc[pd_df['matrix'] == matrix]
    true = df.loc[df['labels'] == 1]
    false = df.loc[df['labels'] == 0]
    sns.distplot(false['scores'],hist=False,label='FPR')
    sns.distplot(true['scores'],hist=False,label='TPR')
    plt.title('Distributions for {0} at {1}%'.format(matrix,threshold))
    plt.legend()

    # save plot to image directory
    outdir = pathlib.Path('images')
    outfile = outdir / 'Dist_{0}_{1}.png'.format(matrix,threshold)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outfile),format='png')
    plt.close()

    # all_labels = label_thresholded_scores(threshold,pos_scores+neg_scores)
def label_thresholded_scores(threshold,scores):
    '''
    if the score is over the threshold label 1
    '''
    return [1 if x > threshold else 0 for x in scores]

# def get_threshold_location(pd_scores,threshold):
# My attempt to graph the distriubtions of false and tru positives
#     negatives = pd_scores.loc[pd_scores['labels']==threshold]
#     subset_index = int(len(negatives['scores'])-(len(negatives['scores'])*(threshold/100))-1)
#     print(subset_index)
#     print(subset_index)
#     print(sorted(negatives['scores']))
#     thresh_loc = sorted(negatives['scores'])[subset_index]
#     print(thresh_loc)
#     return thresh_loc
    # pd_stats_0,pd_scores_0 = get_score_dfs(gap_ext,gap,matrix,pos_pairs[:8],neg_pairs[:8],0)
    # thresh_loc_0 = get_threshold_location(pd_scores_0,0)
    # plot_distributions(pd_scores_0,matrix,0)

def question_two_pt_one(pd_scores,gap_ext,gap,matrix,pos_pairs,neg_pairs,threshold):
    '''
    modify starting score matrix
    '''
    # read in cost matrix
    cost_matrix = read_in_cost_matrix(matrix)

    pos_pairs = [line.strip() for line in open('pairs/Pospairs.txt','r')]
    neg_pairs = [line.strip() for line in open('pairs/Negpairs.txt','r')]

    # num_iterations = 20 # number of iterations over which to optimize
    # num_permutations = 10 # number of permuted matrices to make
    # for i in range(num_iterations):
    #     collect_matrices = []
    #     cost_values = cost_matrix.values
    #     for j in range(num_permutations):
    #         cost_random = np.random.permutation(cost_values)
    #         cost_sym = np.tril(cost_random) + np.tril(cost_random, -1).T
    #         pd_stats_0,pd_scores_0 = get_score_dfs(gap_ext,gap,cost_sym,pos_pairs[:8],neg_pairs[:8],0)
            # evaluate the score
            # if the matrix is a reasonable improvment, add to the collection and shuffle
            # else discard
