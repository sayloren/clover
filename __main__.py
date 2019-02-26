import pandas as pd
import argparse
import os
from alignment import smith_waterman,read_fasta,get_score_dfs
from question import question_one_pt_one,question_one_pt_two,question_one_pt_three,question_two_pt_one

def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument("-t","--threshold",type=int,default="70",help='the percentage at which to threshold true positives')
    parser.add_argument("-g","--gap",type=int,default="8",help='the gap opening size to run as a range from 1') # 21
    parser.add_argument("-e","--gapext",type=int,default="3",help='the gap extention penalty to run as a range from 1') # 6
    parser.add_argument("-m","--matrix",type=str,default="PAM250.txt",help='if runing a single pass, which matrix to use') # 6
    parser.add_argument("-s","--single",action='store_true',help='if want to run for a single set of gap and ext conditions rather than a range') # 15,3
    parser.add_argument("-i","--inmatrix",action='store_true',help='if data file should be read in, other wise the algorithm is run and the matrix generated')
    return parser.parse_args()

def main():
    args = get_args()

    # 1) get pairs
    pos_pairs = [line.strip() for line in open('pairs/Pospairs.txt','r')]
    neg_pairs = [line.strip() for line in open('pairs/Negpairs.txt','r')]

    # if arg given, read in matrix from file
    if args.inmatrix:
        pd_collect = pd.read_csv('rate_run.csv',sep='\t')
        pd_scores = pd.read_csv('score_run.csv',sep='\t')
    else:
        # initialze the collection data frame
        collect = []
        collect_scores = []

        # if just want to run for a single gap and extention value
        if args.single:
            for m in os.listdir('matrices'):
                if m.endswith('.txt'):
                    out,pd_scores = get_score_dfs(args.gapext,args.gap,m,pos_pairs,neg_pairs,args.threshold)
                    collect.append(out)
                    collect_scores.append(pd_scores)

        # if just want to run for a range of gap and extention values
        else:
            for m in os.listdir('matrices'):
                if m.endswith('.txt'):
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
    # best_gap,best_ext = question_one_pt_one(pd_collect)
    # print('for BLOSUM50, best gap is {0}, best ext is {1}'.format(best_gap,best_ext))

    # question 1.2
    # gap/ext values from above, which matrix has best score
    # best_gap,best_ext = 8,3
    # best_matrix = question_one_pt_two(pd_collect,pd_scores,best_gap,best_ext)
    # print('For gap {0} and ext {1} the best performing matrix is {2}'.format(best_gap,best_ext,best_matrix))

    # question 1.3
    # question_one_pt_three(pd_scores,best_matrix)

    # question 2.1
    question_two_pt_one(pd_scores,args.gapext,args.gap,args.matrix,pos_pairs,neg_pairs,args.threshold)

if __name__ == "__main__":
    main()
