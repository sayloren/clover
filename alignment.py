# script with implmenting smith-waterman from
# https://tiefenauer.github.io/blog/smith-waterman/

def make_scoring_matrix(string_a,string_b,match_score,gap_cost):
    '''
    make a matrix of dimensions string_a x string_b
    populated by the match score between the two strings
    with allowances for gaps
    '''

    # initialze a matrix with 0's
    matrix_init = np.zeros((len(string_a) + 1, len(string_b) + 1), np.int)

    # for row and column along matrix
    for i, j in itertools.product(range(1, matrix_init.shape[0]), range(1, matrix_init.shape[1])):

        # match is the previous value + the match score if those values are the same in a and b, - other wise
        match = matrix_init[i - 1, j - 1] + (match_score if string_a[i - 1] == string_b[j - 1] else - match_score)

        # delete is the same, but misaligned along one axis by the gap score
        delete = matrix_init[i - 1, j] - gap_cost

        # insert is the same, but misaligned along the other axis by the gap score
        insert = matrix_init[i, j - 1] - gap_cost

        # the cell becomes the max value of these options
        matrix_init[i, j] = max(match, delete, insert, 0)

    return matrix_init

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

def smith_waterman(string_a, string_b, match_score=3, gap_cost=2):
    '''
    perform smith waterman algorithm, getting the traceback along the
    matrix of scored values between strings a and b
    '''

    # make all the strings upper case letters
    string_a, string_b = string_a.upper(), string_b.upper()

    # create the scoring matrix for strings a and b
    score_matrix = make_scoring_matrix(string_a, string_b, match_score, gap_cost)

    # get the optimal trace back through the scoring matrix and find the starting point in the string
    string_b_, position = optimal_traceback(score_matrix, string_b)

    # return the starting point in the string, and the length of the matched string b
    return position, position + len(string_b_)
