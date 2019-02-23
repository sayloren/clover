import pandas as pd
import argparse
import os
from itertools import compress
from alignment import smith_waterman,read_fasta
from question import question_one_pt_one,question_one_pt_two,question_one_pt_three

def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument("-t","--threshold",type=int,default="70",help='the percentage at which to threshold true positives')
    parser.add_argument("-g","--gap",type=int,default="21",help='the gap opening size to run as a range from 1')
    parser.add_argument("-e","--gapext",type=int,default="6",help='the gap extention penalty to run as a range from 1')
    parser.add_argument("-s","--single",action='store_true',help='if want to run for a single set of gap and ext conditions rather than a range') # 15,3
    parser.add_argument("-i","--inmatrix",action='store_true',help='if data file should be read in, other wise the algorithm is run and the matrix generated')
    return parser.parse_args()

def get_score_dfs(gap_ext,gap,m,pos_pairs,neg_pairs,threshold):
    # 2) get matrix scoring fasta pair for each matrix
    pos_out = [smith_waterman(x.split()[0], x.split()[1],m,gap,gap_ext) for x in pos_pairs]
    pos_scores,pos_min = [p[0] for p in pos_out],[p[1] for p in pos_out]
    neg_out = [smith_waterman(x.split()[0], x.split()[1],m,gap,gap_ext) for x in neg_pairs]
    neg_scores,neg_min = [n[0] for n in neg_out],[n[1] for n in neg_out]

    # the index of the value from which to subset the sorted true positive list
    subset_index = int(len(pos_scores)-(len(pos_scores)*threshold/100)-1)
    threshold = sorted(pos_scores)[subset_index]
    norm_threshold = sorted(pos_min)[subset_index]

    # 3) get the false and true positive rates
    false_pos = sum([x > threshold for x in neg_scores])/len(neg_scores)
    true_pos = sum([x > threshold for x in pos_scores])/len(pos_scores)

    false_norm = sum(list(compress(neg_min, [x > norm_threshold for x in neg_min])))/len(neg_min)
    true_norm = sum(list(compress(pos_min, [x > norm_threshold for x in pos_min])))/len(pos_min)

    # get the gap, gap extention, false and true positives
    out = [true_pos,false_pos,gap,gap_ext,m,true_norm,false_norm]

    # collect the scores and other values for making roc curve
    all_score = pos_scores + neg_scores
    all_norm = pos_min + neg_min
    all_matrix = [m]*len(pos_scores) + [m]*len(neg_scores)
    all_labels = [1]*len(pos_scores) + [0]*len(neg_scores)
    all_gaps = [gap]*len(pos_scores) + [gap]*len(neg_scores)
    all_ext = [gap_ext]*len(pos_scores) + [gap_ext]*len(neg_scores)
    pd_scores = pd.DataFrame({'scores':all_score,'norm_scores':all_norm,'matrix':all_matrix,'labels':all_labels,'gap':all_gaps,'ext':all_ext})

    return out,pd_scores

def main():
    args = get_args()

    # 1) get pairs
    pos_pairs = [line.strip() for line in open('pairs/Pospairs.txt','r')]
    neg_pairs = [line.strip() for line in open('pairs/Negpairs.txt','r')]

    # initialze the collection data frame
    collect = []
    collect_scores = []

    # if arg given, read in matrix from file
    if args.inmatrix:
        pd_collect = pd.read_csv('rate_run.csv',sep='\t')
        pd_scores = pd.read_csv('score_run.csv',sep='\t')
    else:
        # iterate through matrices, gap and gap extention sizes
        for m in os.listdir('matrices'):

            # if just want to run for a single gap and extention value
            if args.single:
                out,pd_scores = get_score_dfs(args.gapext,args.gap,m,pos_pairs,neg_pairs,args.threshold)
                collect.append(out)
                collect_scores.append(pd_scores)

            # if just want to run for a range of gap and extention values
            else:
                for gap_ext in range(1,args.gapext):
                    for gap in range(1,args.gap):
                        out,pd_scores = get_score_dfs(gap_ext,gap,m,pos_pairs,neg_pairs,args.threshold)
                        collect.append(out)
                        collect_scores.append(pd_scores)

        pd_collect = pd.DataFrame(collect,columns=['true','false','gap','ext','matrix','t_min','f_min'])
        pd_collect.to_csv('rate_run.csv',sep='\t',index=False)

        pd_scores = pd.concat(collect_scores,axis=0)
        pd_scores.to_csv('score_run.csv',sep='\t',index=False)

    # question 1.1
    # get gap/ext values for best false positive rate in BLOSUM50 matrices
    best_gap,best_ext = question_one_pt_one(pd_collect)
    print('for BLOSUM50, best gap is {0}, best ext is {1}'.format(best_gap,best_ext))

    # question 1.2
    # gap/ext values from above, which matrix has best score
    best_matrix = question_one_pt_two(pd_collect,pd_scores,best_gap,best_ext)
    print('For gap {0} and ext {1} the best performing matrix is {2}'.format(best_gap,best_ext,best_matrix))

    # question 1.3
    question_one_pt_three(pd_scores,best_matrix)

if __name__ == "__main__":
    main()
