"""Calculate distances between protein sequences based on the W-metric (Wm).

Reference:
    1. Vinga, Gouveia-Oliveira, Almeida. (2004) Bioinformatics. 20(2):206-215
       doi: 10.1093/bioinformatics/btg392

W-metric includes one-tuple composition information (the difference
in amino acid frequencies between two proteins) and weights from
the scoring matrices used in alignment methods.

"""
import itertools
import numpy as np


def count_seq_chars(seq, alphabet):
    """Count characters from given alphabet that are present in sequence.

    Args:
       seq (str): sequence
       alphabet (str/list): list of allowed characters

    Returns:
       A list of characters' counting occurrences.

    Examples:
       >>> alphabet = 'ACDEFGHIKLMNPQRSTVWY'
       >>> seq = 'MKSTGWHFSG'
       >>> print(count_seq_chars(seq, alphabet))
       [0, 0, 0, 0, 1, 2, 1, 0, 1, 0, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0]

    """
    l = [0 for c in alphabet]
    for i, c in enumerate(alphabet):
        l[i] += seq.count(c)
    return l


def freq_seq_chars(counts):
    """Calculate frequencies of characters (symbols) in a sequence based on
    characters' counts.

    Args:
       counts (list): result of the `count_seq_chars` function
       seqlen (int):  length of a sequence

    Returns:
       A list of frequencies corresponding to alphabet

    Examples:
        >>> l = [0, 0, 0, 0, 1, 2, 1, 0, 1, 0, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0]
        >>> print(freq_seq_chars(l))
        [0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.1,
         0.0, 0.1, 0.0, 0.1, 0.0, 0.0, 0.0,
         0.0, 0.2, 0.1, 0.0, 0.1, 0.0]

    """
    seqlen = float(sum(counts))
    return [c / seqlen for c in counts]


def freq_seqs_chars(seq_records, alphabet):
    """Calculate frequencies of characters from given alphabet
    for multiple sequences (stored as seq_records object).

    Args:
       seq_records (obj): instance of SeqRecords()
       alphabet (list): list of allowed characters

    Returns:
       numpy.ndarray
    """
    l = []
    for i in range(seq_records.count):
        seq = seq_records.seq_list[i]
        counts = count_seq_chars(seq, alphabet)
        freq = freq_seq_chars(counts)
        l.append(freq)
    return np.array(l)


class Distance:
    """Combine vector with a distance function.

    Attributes:
        freqs (ndarray): matrix of sequence-representing vectors
        matrix (ndarray): substitution matrix for amino acid changes

    """

    def __init__(self, seq_records, matrix):
        """Create a instance of Distance.

        Args:
            seq_records (obj: seqrecords.SeqRecords)
            matrix (obj: utils.data.subsmat.SubsMat)

        Examples:
            >>> from .utils.data import subsmat
            >>> from .utils.seqrecords import SeqRecords
            >>> matrix = subsmat.get('blosum62')
            >>> seq_records = SeqRecords()
            >>> seq_records.add('seq1', 'MKSTGWHF')
            >>> seq_records.add('seq2', 'MKSSSSTGWGWG')
            >>> seq_records.add('seq3', 'MKSTLKNGTEQ')

            >>> dist = Distance(seq_records, matrix)

        """

        self.freqs = freq_seqs_chars(seq_records, matrix.alphabet_list)
        self.matrix = matrix

    def pairwise_distance(self, seqnum1, seqnum2):
        """Compute W-metric between two proteins.

        The distance is defined by one-tuple frequencies
        fx and fy of two proteins, weighted by matrix W.

        Formula:
        d^{w} = \sum_{i\in A}\sum_{j\in A}(f_{i}^{X}-f_{i}^{y})
        \cdot (f_{j}^{X}-f_{j}^{y})\cdot w_{ij}

        """
        weights = self.matrix.data
        freqs1 = self.freqs[seqnum1]
        freqs2 = self.freqs[seqnum2]
        f = freqs1 - freqs2
        m = np.outer(f, f) * self.matrix.data
        return np.sum(m)


def main():
    from .utils import distmatrix
    from .utils.data import subsmat
    from .utils.seqrecords import SeqRecords

    matrix = subsmat.get('blosum62')

    seq_records = SeqRecords()
    seq_records.add('seq1', 'MKSTGWHF')
    seq_records.add('seq2', 'MKSSSSTGWGWG')
    seq_records.add('seq3', 'MKSTLKNGTEQ')

    dist = Distance(seq_records, matrix)

    # print dist.pairwise_distance(0, 1)
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()


if __name__ == '__main__':
    main()
