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

def plot_roc_curve(pd_collect):
    '''
    plot receiver operating characteristic (ROC) curve for each matrix
    '''
    # each positive pair value make 1
    # each negative pair value make 0
    # y_true = [1]*len(pos_pairs) + [0]*len(neg_pairs)
    # make both into list

    # fpr,tpr,threshold = [roc_curve(bin_list,matrix) for x in pd_collect]
    # area_under_roc_curve = [auc(fpr,tpr) for x in pd_collect]
    # plt.plot(fpr,tpr,label='matrix, auc = {0}'.format(area_under_roc_curve))
    # plt.title('ROC curves')
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.xlabel('False Positive Rate')
    # plt.ylabel('True Positive Rate')

def question_one_pt_two(pd_collect,best_gap,best_ext):
    '''
    gap/ext values from question one pt one, which matrix has best false positive score
    '''
    # get a subset for those with the gap/ext values from the best false positive values from BLOSUM50
    best_gap_and_ext = pd_collect.loc[(pd_collect['gap']==best_gap) & (pd_collect['ext']==best_ext)]

    # get the matrix that has the best false positive rates
    max_index = best_gap_and_ext['false'].idxmax()

    # get the matrix that from that index, 4 for matrix column
    best_matrix = best_gap_and_ext['matrix'][max_index]

    # plot the roc curves for each matrix
    plot_roc_curve(pd_collect)

    return best_matrix

def question_one_pt_three(pd_collect):
    '''
    normalize smith and waterman score by length of shortest pair sequence
    roc curves for best matrix and normalized scores
    '''
    plot_roc_curve(pd_collect)
