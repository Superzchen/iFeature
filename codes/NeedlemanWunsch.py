#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys
from collections import namedtuple

MatrixCell = namedtuple('MatrixCell', 'score h v')


def print_matrix(m, seq):
    """
    Printing sequences and matrix
    :param m: [][]int
    :param seq: [2]string
    :return:
    """
    print('{:10}'.format("")),
    for char in seq[1]:
        print('{:4}'.format(char)),
    print
    for i in range(len(m)):
        if not i:
            print('{:2}'.format("")),
        else:
            print('{:2}'.format(seq[0][i - 1])),
        for j in range(len(m[i])):
            print('{:4}'.format(m[i][j].score)),
        print

def print_matrix_h(m, seq):
    """
    Printing sequences and matrix
    :param m: [][]int
    :param seq: [2]string
    :return:
    """
    print('{:10}'.format("")),
    for char in seq[1]:
        print('{:4}'.format(char)),
    print
    for i in range(len(m)):
        if not i:
            print('{:2}'.format("")),
        else:
            print('{:2}'.format(seq[0][i - 1])),
        for j in range(len(m[i])):
            print('{:4}'.format(m[i][j].h)),
        print

def print_matrix_v(m, seq):
    """
    Printing sequences and matrix
    :param m: [][]int
    :param seq: [2]string
    :return:
    """
    print('{:10}'.format("")),
    for char in seq[1]:
        print('{:4}'.format(char)),
    print
    for i in range(len(m)):
        if not i:
            print('{:2}'.format("")),
        else:
            print('{:2}'.format(seq[0][i - 1])),
        for j in range(len(m[i])):
            print('{:4}'.format(m[i][j].v)),
        print


def gap_line(gap, length):
    if length > 0:
        return gap[0] + (length - 1) * gap[1]
    else:
        return 0


def matrix_filling_NW(seq, s_matrix, gap):
    """
    Filling matrix according to Needleman-Wunsch algorithm
    :param seq: [2]string
    :param s_matrix: dict( char -> dict( char -> int))
    :param gap: [2]int
    :return:  [2]string, int
    """
    f_matrix = []  # Crating F-matrix, init with 0s
    neg_inf = float('-inf')
    for x in range(len(seq[0]) + 1):
        f_matrix.append([MatrixCell(0, neg_inf, neg_inf) for x in range(len(seq[1]) + 1)])

    for i in range(1, len(f_matrix)):  # Filling first row and column
        f_matrix[i][0] = MatrixCell(gap_line(gap, i), neg_inf, gap_line(gap, i))
    for j in range(1, len(f_matrix[0])):
        f_matrix[0][j] = MatrixCell(gap_line(gap, j), gap_line(gap, j), neg_inf)

    for i in range(1, len(f_matrix)):  # Filling all other cells
        #sys.stdout.write('-')
        for j in range(1, len(f_matrix[i])):
            d = f_matrix[i - 1][j - 1].score + s_matrix[seq[0][i - 1]][seq[1][j - 1]]
            h = max(f_matrix[i][j - 1].score + gap[0], f_matrix[i][j - 1].h + gap[1])
            v = max(f_matrix[i - 1][j].score + gap[0], f_matrix[i - 1][j].v + gap[1])
            score = max(d, h, v)  # Choosing max of variants
            f_matrix[i][j] = MatrixCell(score, h, v)
    #print('\nF-Matrix is built.')
    return result_seq(f_matrix, seq, s_matrix, gap)


def result_seq(f_matrix, seq, s_matrix, gap):
    """
    Finding trace in F matrix and getting resulting sequences
    :param f_matrix: [][]MatrixCell
    :param seq: [2]string
    :param s_matrix: dict( char -> dict( char -> int))
    :param gap: [2]int
    :return: [2]string, int
    """
    res1 = ""
    res2 = ""
    i = len(seq[0])
    j = len(seq[1])
    while i > 0 or j > 0:
        #sys.stdout.write('-')
        if i > 0 and j > 0 and f_matrix[i][j].score == f_matrix[i - 1][j - 1].score + s_matrix[seq[0][i - 1]][seq[1][j - 1]]:
            i -= 1
            j -= 1
            res1 = seq[0][i] + res1
            res2 = seq[1][j] + res2
        elif i > 0 and f_matrix[i][j].score == f_matrix[i][j].v:
            i -= 1
            res1 = seq[0][i] + res1
            res2 = "-" + res2
        elif j > 0 and f_matrix[i][j].score == f_matrix[i][j].h:
            j -= 1
            res1 = "-" + res1
            res2 = seq[1][j] + res2
    #print('\nSequences are found.')
    return res1, res2, f_matrix.pop().pop().score
