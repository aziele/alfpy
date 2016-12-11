"""Distance methods between two boolean vectors (representing word
occurrences).

References:
    1. SciPy, https://www.scipy.org

"""

import numpy as np

from .utils import distance


def _nbool_correspond_ft_tf(u, v):
    """Function used by some distance methods (in Distance class).
    Based on: https://github.com/scipy/scipy

    Args:
        u (numpy.ndarray) : boolean vector, shape: (N, 1)
        v (numpy.ndarray) : as above

    Returns:
        tuple of two numbers

    Examples:
        >>> u = np.array([True, False, True])
        >>> v = np.array([True, True, False])
        >>> print(_nbool_correspond_ft_tf(u, v))
        (1, 1)

    """
    not_u = ~u
    not_v = ~v
    nft = (not_u & v).sum()
    ntf = (u & not_v).sum()
    return (nft, ntf)


def _nbool_correspond_all(u, v):
    """Function used by some distance methods (in Distance class).
    Based on: https://github.com/scipy/scipy

    Args:
        u (numpy.ndarray) : bool, shape: (N, )
        v (numpy.ndarray) : as above

    Returns:
        tuple of four numbers

    Examples:
        >>> u = np.array([True, False, True])
        >>> v = np.array([True, True, False])
        >>> print(_nbool_correspond_all(u, v))
        (0, 1, 1, 1)

    """
    not_u = ~u
    not_v = ~v
    nff = (not_u & not_v).sum()
    nft = (not_u & v).sum()
    ntf = (u & not_v).sum()
    ntt = (u & v).sum()
    return (nff, nft, ntf, ntt)


class Distance(distance.Distance):
    """Combine vector boolean data (numpy.ndarray) with distance method.

    """

    def pwdist_dice(self, seq1idx, seq2idx):
        """Compute the Dice dissimilarity (Sorensen-Dice coefficient)
        between two boolean 1-D arrays.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]
        ntt = (u & v).sum()
        (nft, ntf) = _nbool_correspond_ft_tf(u, v)
        return float(ntf + nft) / float(2.0 * ntt + ntf + nft)

    def pwdist_yule(self, seq1idx, seq2idx):
        """Compute the Yule dissimilarity between two boolean 1-D arrays.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]
        (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)
        return float(2.0 * ntf * nft) / float(ntt * nff + ntf * nft)

    def pwdist_rogerstanimoto(self, seq1idx, seq2idx):
        """Compute the Rogers-Tanimoto dissimilarity between two boolean
        1-D arrays.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]
        (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)
        r = float(2.0 * (ntf + nft)) / float(ntt + nff + (2.0 * (ntf + nft)))
        return r

    def pwdist_russellrao(self, seq1idx, seq2idx):
        """Compute the Russell-Rao dissimilarity between two boolean 1-D arrays.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]

        ntt = (u & v).sum()
        return float(len(u) - ntt) / float(len(u))

    def pwdist_sokalmichener(self, seq1idx, seq2idx):
        """Compute the Sokal-Michener dissimilarity
        between two boolean 1-D arrays.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]
        ntt = (u & v).sum()
        nff = (~u & ~v).sum()
        (nft, ntf) = _nbool_correspond_ft_tf(u, v)
        return float(2.0 * (ntf + nft)) / float(ntt + nff + 2.0 * (ntf + nft))

    def pwdist_sokalsneath(self, seq1idx, seq2idx):
        """Compute the Sokal-Sneath dissimilarity
        between two boolean 1-D arrays.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]
        ntt = (u & v).sum()

        (nft, ntf) = _nbool_correspond_ft_tf(u, v)
        denom = ntt + 2.0 * (ntf + nft)
        if denom == 0:
            raise ValueError('Sokal-Sneath dissimilarity is not defined for '
                             'vectors that are entirely false.')
        return float(2.0 * (ntf + nft)) / denom

    def pwdist_jaccard(self, seq1idx, seq2idx):
        """Compute the Jaccard-Needham dissimilarity
        between two boolean 1-D arrays.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]
        dist = (np.double(np.bitwise_and((u != v),
                np.bitwise_or(u != 0, v != 0)).sum()) /
                np.double(np.bitwise_or(u != 0, v != 0).sum()))
        return dist

    def pwdist_hamming(self, seq1idx, seq2idx):
        """Compute the Hamming distance between two 1-D arrays.

        The Hamming distance between 1-D arrays `u` and `v`, is simply the
        proportion of disagreeing components in `u` and `v`.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]
        return (u != v).mean()

    def pwdist_kulsinski(self, seq1idx, seq2idx):
        """Compute the Kulsinski dissimilarity between two boolean 1-D arrays.

        Returns:
            distance value (double)

        """
        u = self[seq1idx]
        v = self[seq2idx]
        n = float(len(u))
        (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v)
        return (ntf + nft - ntt + n) / (ntf + nft + n)


def main():
    from .utils.seqrecords import SeqRecords
    from . import word_vector
    from . import word_pattern
    from .utils import distmatrix

    seq_records = SeqRecords()
    seq_records.add('seq1', 'MKSTGWHF')
    seq_records.add('seq2', 'MKSSSSTGWGWG')
    seq_records.add('seq3', 'MKSTLKNGTEQ')

    p = word_pattern.create(seq_records.seq_list, 2)
    bools = word_vector.Bools(seq_records.length_list, p)
    dist = Distance(bools, 'jaccard')
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()


if __name__ == '__main__':
    main()
