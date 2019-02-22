import pandas as pd
import argparse
import os
from alignment import smith_waterman,read_fasta
from question import question_one_pt_one,question_one_pt_two

def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument("-t","--threshold",type=int,default="70",help='the percentage at which to threshold true positives')
    parser.add_argument("-g","--gap",type=int,default="21",help='the gap opening size')
    parser.add_argument("-e","--gapext",type=int,default="6",help='the gap extention penalty')
    parser.add_argument("-i","--inmatrix",action='store_true',help='if data file should be read in, other wise the algorithm is run and the matrix generated')
    return parser.parse_args()

def main():
    args = get_args()

    # 1) get pairs
    pos_pairs = [line.strip() for line in open('pairs/Pospairs.txt','r')]
    neg_pairs = [line.strip() for line in open('pairs/Negpairs.txt','r')]

    # initialze the collection data frame
    collect = []

    # if arg given, read in matrix from file
    if args.inmatrix:
        pd_collect = pd.read_csv('paired_runs.csv',sep='\t')
    else:
        # iterate through matrices, gap and gap extention sizes
        for m in os.listdir('matrices'):
            for gap_ext in range(1,args.gapext):
                for gap in range(1,args.gap):

                    # 2) get matrix scoring fasta pair for each matrix
                    pos_scores,pos_min = [smith_waterman(x.split()[0], x.split()[1],m,gap,gap_ext) for x in pos_pairs]
                    neg_scores,neg_min = [smith_waterman(x.split()[0], x.split()[1],m,gap,gap_ext) for x in neg_pairs]

                    # the index of the value from which to subset the sorted true positive list
                    subset_index = int(len(pos_scores)-(len(pos_scores)*args.threshold/100)-1)
                    threshold = sorted(pos_scores)[subset_index]

                    # 3) get the false and true positive rates
                    false_pos = sum([x > threshold for x in neg_scores])/len(neg_scores)
                    true_pos = sum([x > threshold for x in pos_scores])/len(pos_scores)

                    # get the gap, gap extention, false and true positives
                    out = [true_pos,false_pos,gap,gap_ext,m,pos_min,neg_min]
                    collect.append(out)

        pd_collect = pd.DataFrame(collect,columns=['true','false','gap','ext','matrix','t_min','f_min'])
        pd_collect.to_csv('paired_runs.csv',sep='\t')

    # get gap/ext values for best false positive rate in BLOSUM50 matrices
    best_gap,best_ext = question_one_pt_one(pd_collect)
    print(best_gap,best_ext)

    # gap/ext values from above, which matrix has best score
    best_matrix = question_one_pt_two(pd_collect,best_gap,best_ext)
    print(best_matrix)

    # receiver operating characteristic (ROC) curve
    # false pos on x, true pos on y
    # limit x,y to 0-1
    # pd_collect['overall'] = pd_collect['true']+pd_collect['false']




if __name__ == "__main__":
    main()
