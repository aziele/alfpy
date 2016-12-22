"""This module computes distances between DNA/protein sequences based on the
sequence feature, named Base-Base Correlation (BBC).

References:
    1. Liu, Zhi-Hua, et al. (2007) Bioinformatics and Biomedical Engineering,
       ICBBE. The 1st International Conference on. IEEE, 2007.
       doi: 10.1109/ICBBE.2007.98

    2. Liu Z, Meng J, Sun X. (2008) Biochem Biophys Res Commun. 368(2):223-30.
       doi: 10.1016/j.bbrc.2008.01.070.

Todo:
    * handle sequence symbols not included in molecule's alphabet

"""

import numpy as np

from .utils import distance


def base_base_correlation(seq, k, alphabet):
    """Compute the base base correlation (BBC) vector for a sequence.

    Args:
        seq (str) : sequence
        k (int)   : parameter of the BBC. Intuitively, it represents
                    the maximum distance to observe correlation between bases.
        alphabet (str/list) : List of possible characters. This can be used to
                    avoid autodetection of the alphabet in the case where
                    sequences with missing letters are to be compared.

    Returns:
        numpy.ndarray: shape (1, 16) for DNA and (1, 400) for protein.

    Examples:
        >>> print(base_base_correlation('ATGCATGC', 1, 'ATGC'))
        [[
         -0.12547302 -0.12547302  0.2281059   0.17169665  0.01815213
         -0.12547302 -0.12547302  0.04258163  0.04258163  0.17169665
         -0.12547302 -0.12547302 -0.12547302  0.2281059   0.17169665
         -0.12547302
        ]]

    Note:
        A description of the method can be found here:
        http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4272582

        This implementation is generalized for any sequence type.
    """

    s = seq

    if k > len(s) - 2:
        raise Exception("Sequence too short to compute BBC with "
                        "k={}".format(k))

    if alphabet is None:
        alphabet = set(s)
    else:
        s = "".join([c for c in s if c in alphabet])

    alphabet = sorted(list(alphabet))
    alphabet = dict(zip(alphabet, range(len(alphabet))))
    L = len(alphabet)

    # Compute the base probabilities for every character.
    p = np.zeros(L)
    for c in s:
        p[alphabet[c]] += 1
    p /= np.sum(p)
    p.shape = (1, L)

    bbc = np.zeros((L, L))
    for l in range(1, k + 2):
        # Compute $p_{ij}(l)$ representing the probability of
        # observing the bases i and j separated by l "gaps".
        # Compute it for all 16 combinations of alleles.
        l_dist_correlations = np.zeros((L, L))
        for i in range(len(s) - l):
            nuc1 = alphabet[s[i]]
            nuc2 = alphabet[s[i + l]]
            l_dist_correlations[nuc1][nuc2] += 1
        l_dist_correlations /= np.sum(l_dist_correlations)

        # Compute the D_{ij}(l) which is the deviation from
        # statistical independance.
        # $D_{ij}(l) = p_{ij}(l) - p_i p_j$
        D = l_dist_correlations - np.dot(p.T, p)

        bbc += D + (D ** 2 / 2 * np.dot(p.T ** 2, p ** 2)) + D ** 3

    # Flatten the bbc into a 16 feature vector.
    bbc.shape = (1, L * L)
    return bbc


def create_vectors(seq_records, k=10, alphabet="ATGC"):
    """Create BBC's vectors for multiple sequence records.

    Args:
        seq_records (obj SeqRecords)
    """
    data = np.zeros(shape=(seq_records.count, len(alphabet)**2))
    for seqidx, seq in enumerate(seq_records.seq_list):
        vector = base_base_correlation(seq, k=k, alphabet=alphabet)
        data[seqidx] = vector
    return data


class Distance(distance.Distance):

    def __init__(self, vector, disttype='euclid_norm'):
        super(Distance, self).__init__(vector, disttype)


def main():
    from .utils.seqrecords import main
    from .utils import distmatrix
    seq_records = main()
    vector = create_vectors(seq_records, 10, alphabet="ATGC")
    dist = Distance(vector)
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()


if __name__ == '__main__':
    main()
