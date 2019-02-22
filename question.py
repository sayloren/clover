import pandas as pd
from sklearn.metrics import roc_curve, auc

def question_one_pt_one(pd_collect):
    '''
    get which BLOSUM50 gap/ext combo run has the best false positive rate
    '''
    # get all the BLOSUM50 matrices
    blosum_fifty = pd_collect.loc[pd_collect['matrix'] == 'BLOSUM50']

    # get the index of the false positive with the best score
    max_index = blosum_fifty['false'].argmax()

    # get the gap and extension values for that best false positive score
    # 2 for gap column, 3 for ext column
    best_gap,best_ext = blosum_fifty.iloc[2:3,max_index]
    return best_gap,best_ext

def plot_roc_curve(pd_collect):
    '''
    plot receiver operating characteristic (ROC) curve for each matrix
    '''


# y_true = [1]*len(pos_pairs) + [0]*len(neg_pairs)
# fpr = dict()
# tpr = dict()
# roc_auc = dict()
#
# for path in scoring_files:
#     fpr[path], tpr[path], _ = roc_curve(y_true, scoring_results[path])
#     roc_auc[path] = auc(fpr[path], tpr[path])
#
#
# colors = ['deeppink','aqua', 'darkorange', 'cornflowerblue']
# for path, color in zip(scoring_files, colors):
#     plt.plot(fpr[path], tpr[path], color=color, lw=2,
#              label='ROC curve of class {0} (area = {1:0.2f})'
#              ''.format(path, roc_auc[path]))
#
# plt.plot([0, 1], [0, 1], 'k--', lw=2)
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver operating characteristic to all scoring matrices')
# plt.legend(loc="lower right")
# plt.show()




def question_one_pt_two(pd_collect,best_gap,best_ext):
    '''
    gap/ext values from question one pt one, which matrix has best false positive score
    '''
    # get a subset for those with the gap/ext values from the best false positive values from BLOSUM50
    best_gap_and_ext = pd.collect.loc[(pd_collect['gap']==best_gap) & (pd_collect['ext']==best_ext)]

    # get the matrix that has the best false positive rates
    max_index = best_gap_and_ext['false'].argmax()

    # get the matrix that from that index, 4 for matrix column
    best_matrix = blosum_fifty.iloc[4,max_index]

    # plot the roc curves for each matrix
    plot_roc_curve(pd_collect)

    return best_matrix

def question_one_pt_three():
    '''
    normalize smith and waterman score by length of shortest pair sequence
    roc curves for best matrix and normalized scores
    '''
