import pandas as pd
from sklearn.metrics import roc_curve, auc

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

def plot_roc_curve(scores,labels,name,file):
    '''
    plot receiver operating characteristic (ROC) curve for each matrix
    '''
    fpr,tpr,threshold = roc_curve(labels,scores)
    area_under_roc_curve = auc(fpr,tpr)
    plt.plot(fpr,tpr,label='auc = {1}'.format(area_under_roc_curve))
    plt.title('ROC curves for {0}'.format(name))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

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

    # get the best matrix from the scores matrix
    best_scores = pd_scores.loc[pd_scores['matrix'] == best_matrix]

    # get the scores and labels
    labels = best_scores['scores']
    scores = best_scores['labels']

    plot_roc_curve(scores,labels,'best performing matrix, {0}'.format(best_matrix),'Best_'.format(best_matrix))

    return best_matrix

def question_one_pt_three(pd_scores,best_matrix):
    '''
    normalize smith and waterman score by length of shortest pair sequence
    roc curves for best matrix, normalized scores
    '''
    # get the best matrix from the scores matrix
    best_scores = pd_scores.loc[pd_scores['matrix'] == best_matrix]

    # get the scores and labels
    labels = best_scores['norm_scores']
    scores = best_scores['labels']

    plot_roc_curve(scores,labels,'normalized scores, {0}'.format(best_matrix),'Norm_'.format(best_matrix))
