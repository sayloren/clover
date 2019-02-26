import pandas as pd
import argparse
import os
from itertools import compress
from alignment import smith_waterman,read_fasta
from question import question_one_pt_one,question_one_pt_two,question_one_pt_three

def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument("-t","--threshold",type=int,default="70",help='the percentage at which to threshold true positives')
    parser.add_argument("-g","--gap",type=int,default="8",help='the gap opening size to run as a range from 1') # 21
    parser.add_argument("-e","--gapext",type=int,default="3",help='the gap extention penalty to run as a range from 1') # 6
    parser.add_argument("-m","--matrix",type=str,default="PAM250",help='if runing a single pass, which matrix to use') # 6
    parser.add_argument("-s","--single",action='store_true',help='if want to run for a single set of gap and ext conditions rather than a range') # 15,3
    parser.add_argument("-i","--inmatrix",action='store_true',help='if data file should be read in, other wise the algorithm is run and the matrix generated')
    return parser.parse_args()

def run_sw_and_split_out(pairs,m,gap,gap_ext):
    '''
    run the smith waterman alg and split the scoring and source matrices into
    two outputs
    '''
    out = [smith_waterman(x.split()[0], x.split()[1],m,gap,gap_ext) for x in pairs]
    return [x[0] for x in out],[x[1] for x in out],[x[2] for x in out]

def get_rates(threshold,scores):
    '''
    get the false and true positive rates
    '''
    return sum([x > threshold for x in scores])/len(scores)

def get_norm_rates(threshold,min):
    '''
    get the normalized false and true positive rates
    '''
    return sum(list(compress(min, [x > threshold for x in min])))/len(min)

def make_columns(value_a,value_b,length_a,length_b):
    '''
    make columns for output data frame
    '''
    return [value_a]*len(length_a) + [value_b]*len(length_b)

def get_score_dfs(gap_ext,gap,m,pos_pairs,neg_pairs,threshold):
    '''
    make a df from the scoring matrices
    '''
    pos_scores,pos_min,pos_align = run_sw_and_split_out(pos_pairs,m,gap,gap_ext)
    neg_scores,neg_min,neg_align = run_sw_and_split_out(neg_pairs,m,gap,gap_ext)

    # the index of the value from which to subset the sorted true positive list
    subset_index = int(len(pos_scores)-(len(pos_scores)*threshold/100)-1)

    # 3) get the false and true positive rates
    threshold = sorted(pos_scores)[subset_index]
    false_pos = get_rates(threshold,neg_scores)
    true_pos = get_rates(threshold,pos_scores)

    # get the normalized false and true positive rates
    norm_threshold = sorted(pos_min)[subset_index]
    false_norm = get_norm_rates(norm_threshold,neg_min)
    true_norm = get_norm_rates(norm_threshold,pos_min)

    # get the gap, gap extention, false and true positives
    out = [true_pos,false_pos,gap,gap_ext,m,true_norm,false_norm]

    # collect the scores and other values for making roc curve
    all_score = pos_scores + neg_scores
    all_norm = pos_min + neg_min
    all_matrix = make_columns(m,m,pos_scores,neg_scores)
    all_labels = make_columns(1,0,pos_scores,neg_scores)
    all_gaps = make_columns(gap,gap,pos_scores,neg_scores)
    all_ext = make_columns(gap_ext,gap_ext,pos_scores,neg_scores)
    all_align = make_columns(pos_align,neg_align,pos_scores,neg_scores)
    pd_scores = pd.DataFrame({'scores':all_score,'norm_scores':all_norm,
        'matrix':all_matrix,'labels':all_labels,'gap':all_gaps,'ext':all_ext,'alignment':all_align})

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

        # if just want to run for a single gap and extention value
        if args.single:
            m = args.matrix
            out,pd_scores = get_score_dfs(args.gapext,args.gap,m,pos_pairs,neg_pairs,args.threshold)
            collect.append(out)
            collect_scores.append(pd_scores)

        # if just want to run for a range of gap and extention values
        else:
            for m in os.listdir('matrices'):
                for gap_ext in range(1,args.gapext):
                    for gap in range(1,args.gap):
                        out,pd_scores = get_score_dfs(gap_ext,gap,m,pos_pairs,neg_pairs,args.threshold)
                        collect.append(out)
                        collect_scores.append(pd_scores)

        pd_collect = pd.DataFrame(collect,columns=['true','false','gap','ext','matrix','t_min','f_min'])
        pd_collect.to_csv('rate_run.csv',sep='\t',index=False)

        pd_scores = pd.concat(collect_scores,axis=0)
        pd_scores.to_csv('score_run.csv',sep='\t',index=False)

    # # question 1.1
    # # get gap/ext values for best false positive rate in BLOSUM50 matrices
    # best_gap,best_ext = question_one_pt_one(pd_collect)
    # print('for BLOSUM50, best gap is {0}, best ext is {1}'.format(best_gap,best_ext))
    #
    # # question 1.2
    # # gap/ext values from above, which matrix has best score
    # best_matrix = question_one_pt_two(pd_collect,pd_scores,best_gap,best_ext)
    # print('For gap {0} and ext {1} the best performing matrix is {2}'.format(best_gap,best_ext,best_matrix))
    #
    # # question 1.3
    # question_one_pt_three(pd_scores,best_matrix)

    question_two_pt_one(pd_scores,m):

if __name__ == "__main__":
    main()
