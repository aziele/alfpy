"""This module computes distances between DNA sequences based on the Frequency
Chaos Game Representation (FCGR)

References:
    1. Hatje K, Kollmar M (2012) Front Plant Sci 3: 192.
       doi: 10.3389/fpls.2012.00192


Functions for creating DNA-representing vectors were built upon:
    Cheng J, Cao F, Liu Z. (2013) Mol Biol Evol. 2013 30(5):1032-7.
    doi: 10.1093/molbev/mst021.

"""

import numpy as np

from .utils import distance


def fcgr_vector(dnaseq, word_size):
    """Create a FCGR vector representing a DNA sequence.

    Args:
        dnaseq (str/list): dna sequence
        word_size (int): word size (>= 1)

    Returns:
        list (length equals 4^word_size)

    Examples:
        >>> s = 'ATGCTGATGGATG'
        >>> print(fcgr_vector(s, 1))
        [5, 3, 5]

        >>> print(fcgr_vector(s, 2))
        [1, 0, 1, 0, 0, 0, 4, 0, 2, 2, 0, 0, 1, 3, 0]

    """
    ndata = pow(4, word_size)
    genlen = len(dnaseq)
    CGRs = np.zeros((genlen + 1, 2))

    Apoint = np.array((0.0, 1.0))
    Tpoint = np.array((1.0, 1.0))
    Gpoint = np.array((1.0, 0.0))
    Cpoint = np.array((0.0, 0.0))
    CGRs[0, 0] = 0.5
    CGRs[0, 1] = 0.5
    for i in range(0, genlen):
        if dnaseq[i] == 'A':
            CGRs[i + 1] = 0.5 * (CGRs[i] + Apoint)
        if dnaseq[i] == 'T':
            CGRs[i + 1] = 0.5 * (CGRs[i] + Tpoint)
        if dnaseq[i] == 'G':
            CGRs[i + 1] = 0.5 * (CGRs[i] + Gpoint)
        if dnaseq[i] == 'C':
            CGRs[i + 1] = 0.5 * (CGRs[i] + Cpoint)
    temp = 1.0 / pow(2, word_size)

    vectors = np.zeros(shape=(1, ndata))  # numpy
    vectors = [0.0] * ndata  # list

    for point in CGRs:
        xx = int(point[0] / temp)
        yy = int(point[1] / temp)
        if yy == pow(2, word_size):
            yy = pow(2, word_size) - 1
        vectors[yy * pow(2, word_size) + xx] += 1
    vectors.pop(0)
    return vectors


def create_vectors(seq_records, word_size):
    """Create a matrix of FCGR vectors.

    Args:
        seq_records (obj: SeqRecords)
        word_size (int): word size (>= 1)

    Returns:
        numpy.ndarray

    """
    data = np.zeros(shape=(seq_records.count, pow(4, word_size) - 1))
    for seqidx, seq in enumerate(seq_records.seq_list):
        vector = fcgr_vector(seq, word_size)
        data[seqidx] = vector
    return data


class Distance(distance.Distance):

    def __init__(self, vector, disttype='euclid_norm'):
        return super(Distance, self).__init__(vector, disttype)


def main():
    from .utils.seqrecords import main
    from .utils import distmatrix
    seq_records = main()

    vector = create_vectors(seq_records, 1)
    dist = Distance(vector)
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()


if __name__ == '__main__':
    main()
