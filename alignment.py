import re
import os
import numpy as np
import pandas as pd
import math
import itertools
# https://tiefenauer.github.io/blog/smith-waterman/

def make_scoring_matrix(cost_matrix,string_a,string_b,open_gap,extend_gap):
    '''
    make a matrix of dimensions string_a x string_b
    populated by the match score between the two strings
    '''

    # get the lengths of the strings
    len_a = len(string_a) + 1
    len_b = len(string_b) + 1

    # initialze matrices with 0's
    # first is aligned, second is misaligned with gaps in a, third is misaligned with gaps in b
    matrix_init = np.zeros((len_a, len_b))
    matrix_a_init = np.zeros((len_a, len_b))
    matrix_b_init = np.zeros((len_a, len_b))

    # where there is no sequence to align to set to inf
    matrix_init[[i for i in range(1,len_a)],0] = -(math.inf)
    matrix_init[0,[j for j in range(1,len_b)]] =  -(math.inf)
    matrix_a_init[[i for i in range(1,len_a)],0] =  -(math.inf)
    matrix_b_init[0,[j for j in range(1,len_b)]] =  -(math.inf)

    # gaps in first seq for every position in second for alignment to zero characters
    for j in range(1,len_b):
        matrix_a_init[0,j] = j * extend_gap + open_gap
    for i in range(1,len_a):
        matrix_b_init[i,0] = i * extend_gap + open_gap

    # for row and column along matrix
    for i, j in itertools.product(range(1, matrix_init.shape[0]), range(1, matrix_init.shape[1])):

        # terms
        comb_gap = open_gap-extend_gap
        matrix_a_pos = matrix_init[i,j-1]
        matrix_b_pos = matrix_init[i-1,j]
        str_a_pos = string_a[i-1]
        str_b_pos = string_b[j-1]
        matrix_pre = matrix_init[i-1,j-1]
        matrix_a_pre = matrix_a_init[i-1,j-1]
        matrix_b_pre = matrix_b_init[i-1,j-1]

        # the max value from the cost matrix given string a position i and string b position j, checking for misalignments along either axis
        matrix_a_init[i,j] = max(matrix_a_pos-comb_gap,matrix_a_pos-extend_gap,matrix_a_pos-comb_gap)
        matrix_b_init[i,j] = max(matrix_b_pos-comb_gap,matrix_b_pos-comb_gap,matrix_b_pos-comb_gap)
        matrix_init[i,j] = cost_matrix.loc[str_a_pos,str_b_pos]+max(matrix_pre,matrix_a_pre,matrix_b_pre)

    # get the matrix that has the highest value last element
    max_values = matrix_init[-1,-1],matrix_a_init[-1,-1],matrix_b_init[-1,-1]
    max_matrix_index = max_values.index(max(max_values))
    matrices = [matrix_init,matrix_a_init,matrix_b_init]
    max_matrix = matrices[max_matrix_index]

    return max_matrix

def optimal_traceback(score_matrix, string_b, string_b_='',previous_i=0):
    '''
    recurisively for each substring in string b
    find the optimal traceback path within the scoring matrix
    get each last max occurance starting from the bottom of the scoring matrix
    '''

    # invert the matrix in order to get the last value from both axes
    score_matrix_flip = np.flip(np.flip(score_matrix, 0), 1)

    # get the index of the max value in the score matrix
    i_, j_ = np.unravel_index(score_matrix_flip.argmax(), score_matrix_flip.shape)

    # adjust the index values
    i, j = np.subtract(score_matrix.shape, (i_ + 1, j_ + 1))

    # if 0 it means you've gone through the entire matrix
    if score_matrix[i, j] == 0:
        return string_b_, j

    # other wise add to the string value or add gap string
    string_b_ = string_b[j - 1] + '-' + string_b_ if previous_i - i > 1 else string_b[j - 1] + string_b_

    # recursivley perform through axis b string
    return optimal_traceback(score_matrix[0:i, 0:j], string_b, string_b_, i)

def read_fasta(file):
    '''
    read in the fasta sequence and return sequence and protein name
    '''
    # get protein name
    name = re.split('-|.',file)[1]

    # open file and read all lines, remove first for header
    with open(file,'r') as f:
        lines = f.read().splitlines()
        lines = lines[1:]

    # collect sequences
    sequence = []
    for line in lines:
        sequence = sequence + list(line)

    return sequence

def smith_waterman(file_a,file_b,matrix_file,open_gap,extend_gap):
    '''
    perform smith waterman algorithm, getting the traceback along the
    matrix of scored values between strings a and b
    '''

    # get the sequence for files a and b
    string_a,string_b = read_fasta(file_a),read_fasta(file_b)

    # get the rows that start with comment in the cost matrix
    exclude = [i for i, line in enumerate(open(os.path.join('matrices',matrix_file),"r")) if line.startswith('#')]

    # open matrix file
    m_file = open(os.path.join('matrices', matrix_file), "r")

    # read in the matrix, skipping the rows that are comments
    cost_matrix =  pd.read_fwf(m_file,skiprows = exclude[0:])

    # set the index to be the same as the column names
    cost_matrix.set_index(cost_matrix.columns.values,inplace=True)

    # make sure strings are upper case
    string_a = [x.upper() for x in string_a]
    string_b = [x.upper() for x in string_b]

    # create the scoring matrix for strings a and b
    score_matrix = make_scoring_matrix(cost_matrix,string_a,string_b,open_gap,extend_gap)

    # get the optimal trace back through the scoring matrix and find the starting point in the string
    # string_b_, position = optimal_traceback(score_matrix, string_b)

    # return the starting point in the string, and the length of the matched string b
    #position, position + len(string_b_)

    # get the score of the matrix
    return score_matrix[-1,-1]
