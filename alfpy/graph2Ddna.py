"""
Computing distances between DNA sequences based on graphical representation
of nucleotide sequences.

References:
    1.  Huang G, Zhou H, Li Y, Xu L. (2011) J Theor Biol. 281(1):107-12.
        doi: 10.1016/j.jtbi.2011.04.003.
    2.  Deng M, Yu C, Liang Q et al. (2011) PLoS One. 6(3):e17293.
        doi: 10.1371/journal.pone.0017293.
    3.  Yu C, Liang Q, Yin C et al. (2010) DNA Res. 17(3):155-68.
        doi: 10.1093/dnares/dsq008.

Functions for creating DNA-representing vectors were built upon:
    Cheng J, Cao F, Liu Z. (2013) Mol Biol Evol. 2013 30(5):1032-7.
    doi: 10.1093/molbev/mst021.

"""

import math
import numpy as np

from .utils import distance


def _2DSGraphVector(seq):
    """Create 10-dimensional statistical vector to characterize a DNA sequence.

    Args:
        seq (str/list): DNA sequence

    Returns:
        numpy.ndarray (10,)

    Examples:
        >>> s = 'ATGCCTGACTGNATGAGAGAC'
        >>> print(_2DSGraphVector(s))
        [  2.31 -11.83  12.    -4.89   7.75   0.74  11.5   -3.85  9.5  -2.35]

    """
    nucleotides = [
        ("A", -(3**0.3333333)),
        ("T", 2**0.3333333),
        ("G", -(5**0.5)),
        ("C", 3**0.5)
    ]

    points = np.zeros((len(seq), 2))
    temppoint = np.array((0, 0))

    d = {}
    for nt, val in nucleotides:
        d[nt] = [np.array([1, val]), 0, 0, 0]

    i = 0
    for nt in seq:
        # Include only nucleotides
        if nt in d:
            points[i] = d[nt][0] + temppoint
            d[nt][1] += points[i, 0]
            d[nt][2] += points[i, 1]
            d[nt][3] += 1
            temppoint = points[i]
            i += 1

    peak = points.max(axis=0)[1]
    lowest = points.min(axis=0)[1]
    l = [peak, lowest]

    for nt, _ in nucleotides:
        for i in range(1, 3):
            l.append(d[nt][i] / d[nt][3])

    return np.array(l)


def _2DMGraphVector(seq, n):
    """Create n-dimensional moment vector to characterize a DNA sequence.

    Args:
        seq (str/list): DNA sequence
        n (int):

    Returns:
        numpy.ndarray (n,)

    Examples:
        >>> s = 'ATGCCTGACTGNATGAGAGAC'
        >>> print(_2DMGraphVector(s, 10))
        [21. 13.44  13.44 16.15  21.16   29.01  40.87  58.58  84.99  124.39]

        >>> print(_2DMGraphVector(s, 5))
        [ 21.    13.44  13.44  16.15  21.16]

    """
    nucleotides = [
        ("A", -(3**0.3333333)),
        ("T", 2**0.3333333),
        ("G", -(5**0.5)),
        ("C", 3**0.5)
    ]

    points = np.zeros((len(seq), 2))
    temppoint = np.array((0, 0))

    d = {}
    for nt, val in nucleotides:
        d[nt] = np.array([1, val])

    i = 0
    for nt in seq:
        # Include only nucleotides
        if nt in d:
            points[i] = d[nt] + temppoint
            temppoint = points[i]
            i += 1
    seqlen = i
    l = []
    for k in range(0, n):
        v = np.sum((points[:, 0] - points[:, 1])**k)
        l.append(v / pow(seqlen, k))
    return np.array(l)


def _2DNGraphVector(genomeseq):
    """Create 48-dimensional natural vector to characterize a DNA sequence.

    Args:
        seq (str/list): DNA sequence

    Returns:
        numpy.ndarray (48,)

    Examples:
        >>> s = 'ATGCCTGACTGNATGAGAGAC'
        >>> print(_2DNGraphVector(s), 10)
        [  6.00e+00   4.00e+00   6.00e+00   4.00e+00   1.17e+01   7.00e+00
           1.10e+01   8.75e+00   1.99e+00   9.52e-01   1.51e+00   2.18e+00
          -5.52e-01  -1.76e-01  -4.70e-01  -3.18e-01   5.73e-02   1.66e-02
           4.52e-02   4.10e-02  -5.12e-03  -1.35e-03  -3.82e-03  -4.10e-03
           4.77e-04   1.13e-04   3.35e-04   4.30e-04  -4.41e-05  -9.38e-06
          -2.92e-05  -4.47e-05   4.08e-06   7.81e-07   2.55e-06   4.66e-06
          -3.78e-07  -6.51e-08  -2.23e-07  -4.85e-07   3.50e-08   5.43e-09
           1.94e-08   5.05e-08  -3.24e-09  -4.52e-10  -1.70e-09  -5.26e-09]

    """
    genlen = len(genomeseq)

    nts = 'ATGC'

    idx = {nt: [] for nt in nts}

    for i, nt in enumerate(genomeseq):
        if nt in idx:
            idx[nt].append(i)

    counts = {nt: len(idx[nt]) for nt in nts}
    vector = [counts[nt] for nt in nts]

    means = {nt: np.mean(idx[nt]) for nt in nts}
    for nt in nts:
        vector.append(means[nt])

    templist = {nt: idx[nt] for nt in nts}

    for k in range(2, 12):
        for nt in nts:
            for i in range(0, counts[nt]):
                tempn = pow((idx[nt][i] - means[nt]), k) / \
                    (pow(counts[nt] * genlen, k - 1))
                templist[nt][i] = tempn
            vector.append(np.sum(templist[nt]))
    return np.array(vector)


def create_2DSGraphVectors(seq_records):
    data = np.zeros(shape=(seq_records.count, 10))
    for seqidx, seq in enumerate(seq_records.seq_list):
        vector = _2DSGraphVector(seq)
        data[seqidx] = vector
    return data


def create_2DMGraphVectors(seq_records, n):
    data = np.zeros(shape=(seq_records.count, n))
    for seqidx, seq in enumerate(seq_records.seq_list):
        vector = _2DMGraphVector(seq, n)
        data[seqidx] = vector
    return data


def create_2DNGraphVectors(seq_records):
    data = []
    for seqidx, seq in enumerate(seq_records.seq_list):
        vector = _2DNGraphVector(seq)
        data.append(vector)
    return data


class Distance(distance.Distance):

    def __init__(self, vector, disttype='euclid_norm'):
        return super(Distance, self).__init__(vector, disttype)


def main():
    from .utils.seqrecords import main
    from .utils import distmatrix
    seq_records = main()

    vector = create_2DSGraphVectors(seq_records)
    dist = Distance(vector)
    print(dist.pairwise_distance(2, 2))
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()

    vector = create_2DMGraphVectors(seq_records, 7)
    dist = Distance(vector)
    print(dist.pairwise_distance(2, 2))
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()

    vector = create_2DNGraphVectors(seq_records)
    dist = Distance(vector)
    print(dist.pairwise_distance(2, 2))
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()


if __name__ == '__main__':
    main()
