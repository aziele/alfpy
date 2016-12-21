"""Return Time Distribution distance (RTD)

In contrast to other word-based measures, RTD accounts for the words'
relative orders. Although, originally presented for DNA sequences, the
implemention handles proteins as well.

Return time can be defined as the time required for the reappearance of a
particular state without its appearance within the epoch. The `return time`
in the context of nucleotide sequence can be defined as the number of
nucleotides between the successive appearances of a particular nucleotide(s)
or k-mer. The frequency distribution of those RTs for a particular k-mer is
referred as a return time distribution (RTD) of that k-mer.

References:
    1. Kolekar, Kale, Kulkarni-Kale (2012) Mol Phylogenet Evol 65 510-522
       doi: http://dx.doi.org/10.1016/j.ympev.2012.07.003.

"""

import numpy as np
from .utils import distance


def calc_rtd(word_positions):
    """Compute return time distribution (RTD) of a given word.

    Args:
        word_positions (list) : list of sequence positions of a given word

    Returns:
        mean, stdev (tuple)

    Examples:
        >>> seq = 'CTACACAACTTTGCGGGTAGCCGGAAACATTGTGAATGCGGTGAACA'
        >>> apos = [i for i, nt in enumerate(seq) if nt == 'A']
        >>> print(apos)
        [2, 4, 6, 7, 18, 24, 25, 26, 28, 34, 35, 43, 44, 46]
        >>> print(calc_rtd(apos, 1))
        (3.3846153846153846, 3.1510306381944679)

    """
    l = []
    positions_count = len(word_positions)
    if positions_count < 2:
        return 0.0, 0.0
    for i in range(1, positions_count):
        pos1 = word_positions[i - 1]
        pos2 = word_positions[i]
        pos = pos2 - pos1
        l.append(pos)
    return np.mean(l), np.std(l)


def create_vector(seqcount, pattern):
    """Compute a matrix of sequence-representing RTD vectors

    Args:
        seqcount (int): number of sequences
        pattern (obj: word_pattern.Pattern)

    Returns:
        ndarray: matrix of RTD vectors
                 (shape: number of seqs, doubled number of words)

    """
    words = pattern.pat_list
    data = np.zeros(shape=(seqcount, len(words) * 2))
    for wordidx in range(len(words)):
        for seqidx in pattern.pos_list[wordidx]:
            word_positions = pattern.pos_list[wordidx][seqidx]
            mean, std = calc_rtd(word_positions)
            data[seqidx, wordidx * 2] = mean
            data[seqidx, wordidx * 2 + 1] = std
    return data


class Distance(distance.Distance):
    pass


def main():
    from .utils.seqrecords import main
    from . import word_pattern
    from .utils import distmatrix

    seq_records = main()
    p = word_pattern.create(seq_records.seq_list, 2, True)
    vector = create_vector(seq_records.count, p)
    dist = Distance(vector, 'google')
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()


if __name__ == '__main__':
    main()
