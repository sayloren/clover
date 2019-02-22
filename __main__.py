import pandas as pd
import argparse
from alignment import smith_waterman


def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument("-t","--threshold",type=int,default="70",help='the percentage at which to threshold true positives')
    parser.add_argument("-g","--gap",type=int,default="21",help='the gap opening size')
    parser.add_argument("-e","--gapext",type=int,default="6",help='the gap extention penalty')
    return parser.parse_args()

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
    return best_matrix


def main():
    args = get_args()

    # 1) get pairs
    pos_pairs = [line.strip() for line in open('pairs/Pospairs.txt','r')]
    neg_pairs = [line.strip() for line in open('pairs/Negpairs.txt','r')]

    # initialze the collection data frame
    collect = []

    # iterate through matrices, gap and gap extention sizes
    for m in os.listdir('matrices'):
        for gap_ext in range(1,args.gapext):
            for gap in range(1,args.gap):

                # 2) get matrix scoring fasta pair for each matrix
                pos_scores = [smith_waterman(x.split()[0], x.split()[1],m,gap,gap_ext) for x in pos_pairs]
                neg_scores = [smith_waterman(x.split()[0], x.split()[1],m,gap,gap_ext) for x in neg_pairs]

                # the index of the value from which to subset the sorted true positive list
                subset_index = int(len(pos_scores)-(len(pos_scores)*args.threshold/100)-1)
                threshold = sorted(pos_scores)[subset_index]

                # 3) get the false and true positive rates
                false_pos = sum([x > threshold for x in neg_scores])/len(neg_scores)
                true_pos = sum([x > threshold for x in pos_scores])/len(pos_scores)

                # get the gap, gap extention, false and true positives
                out = [true_pos,false_pos,gap,gap_ext,m]
                collect.append(out)

    pd_collect = pd.DataFrame(collect,columns=['true','false','gap','ext','matrix'])

    # get gap/ext values for best false positive rate in BLOSUM50 matrices
    best_gap,best_ext = question_one_pt_one(pd_collect)

    # gap/ext values from above, which matrix has best score
    best_matrix = question_one_pt_two(pd_collect,best_gap,best_ext)





if __name__ == "__main__":
    main()
