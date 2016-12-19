"""Lempel-Ziv complexity-based distance

Computing distances between DNA/protein sequences based on the
Lempel-Ziv complexity (LZ complexity).

LZ complexity of a finite sequence S is related to the number
of steps required by a production process that builds S.

References:
    1. Lempel A, Ziv. (1976) IEEE Transactions on. 22(1):75-81.
       doi: 10.1109/TIT.1976.1055501

    2. Otu HH, Sayood K. (2003) Bioinformatics. 19(16):2122-30.
       doi: 10.1093/bioinformatics/btg295.

    3. Hohl M, Ragan M. (2007) Systematic Biology. 56(2):206-21.
       doi: 10.1080/10635150701294741

"""

import itertools


def complexity(s):
    """Calculate a measure of algorithmic complexity c
    introduced by Lempel and Ziv.

    As described in:
    Kaspar F, Schuster HG. Phys Rev A. 1987 36(2):842-848.
    doi: http://dx.doi.org/10.1103/PhysRevA.36.842

    Based on:
    http://stackoverflow.com/a/30694008

    Args:
        s (str/list): sequence of any characters

    Returns:
        float
    """
    i, k, l = 0, 1, 1
    k_max = 1
    n = len(s) - 1
    c = 1
    while True:
        if s[i + k - 1] == s[l + k - 1]:
            k = k + 1
            if l + k >= n - 1:
                c = c + 1
                break
        else:
            if k > k_max:
                k_max = k
            i = i + 1
            if i == l:
                c = c + 1
                l = l + k_max
                if l + 1 > n:
                    break
                else:
                    i = 0
                    k = 1
                    k_max = 1
            else:
                k = 1
    return c


def complexity1(seq):
    """Calculate a measure of algorithmic complexity c
    introduced by Lempel and Ziv.

    As described in:
    Otu HH, Sayood K. Bioinformatics. 2003 19(16):2122-30.
    PubMed PMID: 14594718.

    Based on:
    Hohl M, Ragan M. Systematic Biology. 2007. 56(2):206-21

    Args:
        s (str/list): sequence of any characters

    Returns:
        float
    """
    start = 2
    length = len(seq)
    if length < start:
        # Complexity == length.
        return length

    complexity = 2
    hist_len = 1
    pos = start
    last_match = 0

    for index in range(pos, length):
        is_repro, last_match = is_reproducible(
            seq, index, hist_len, last_match)
        if is_repro:
            hist_len += 1
        else:
            hist_len = 1
            complexity += 1

    return complexity


def is_reproducible(seq, index, hist_len, last_match=0):
    """Check whether a string is reproducible.

    A sequence R is reproducible from sequence S (denoted S > R) when
    R can be obtained from S by copying elements from p-th location in
    S to the end of S. For example, AACGT > AACGTCGTCG with p = 3 and
    AACGT > AACGTAC with p = 2.

    Based on:
    Hohl M, Ragan M. Systematic Biology. 2007. 56(2):206-21
    """

    hist_start = index - hist_len
    for i in range(last_match, hist_start):
        # j == hist_len in last iteration
        for j in range(0, hist_len + 1):
            if seq[i + j] != seq[hist_start + j]:
                break
        if j == hist_len:
            # returns index of new last_match
            return True, i

    return False, 0


class Distance:
    """Five measures of sequence distance, introduced by
    Otu and Sayood (2003), based on LZ complexity.
    """

    def __init__(self, seq_records, disttype='d1_star'):
        self.seq_records = seq_records
        # Precomputed L-Z complexity for input sequences
        # as well as all pairwise concatenated sequences.
        self._complexity = self.__precompute_complexity()
        # Set a default distance measure.
        self.set_disttype(disttype)

    def __precompute_complexity(self):
        d = {}
        # Complexity for single input sequences.
        seqs = self.seq_records.seq_list
        for seqidx, seq in enumerate(seqs):
            d[(seqidx,)] = complexity(seq)
        # Complexity for pairwise concatenated sequences.
        for i, j in itertools.combinations(range(self.seq_records.count), 2):
            c1 = d[(i,)]
            c2 = d[(j,)]
            seq12 = seqs[i] + seqs[j]
            c12 = complexity(seq12)
            d[(i, j)] = c12
            seq21 = seqs[j] + seqs[i]
            c21 = complexity(seq21)
            d[(j, i)] = c21
        return d

    def __get_complexity(self, seq1idx, seq2idx):
        # Fetch cached complexity values.
        seqs = self.seq_records.seq_list
        c1 = self._complexity[(seq1idx,)]
        c2 = self._complexity[(seq2idx,)]
        c12 = self._complexity[(seq1idx, seq2idx)]
        c21 = self._complexity[(seq2idx, seq1idx)]
        return c1, c2, c12, c21

    def pwdist_d(self, seq1idx, seq2idx):
        """Given two sequences S and Q, the function d(S,Q) equals:
        d(S,Q) = max\left \{c(SQ)-c(S), c(QS)-c(Q) \right \}
        """
        c1, c2, c12, c21 = self.__get_complexity(seq1idx, seq2idx)
        return max(c12 - c1, c21 - c2)

    def pwdist_d_star(self, seq1idx, seq2idx):
        """Given two sequences S and Q, the function d*(S,Q) equals:
        d^{*}(S,Q) = \frac{max\left \{ c(SQ-c(S), c(QS)-c(Q)
        \right \}}{max\left \{ c(S),c(Q) \right \}}
        """
        c1, c2, c12, c21 = self.__get_complexity(seq1idx, seq2idx)
        return float(max(c12 - c1, c21 - c2)) / max(c1, c2)

    def pwdist_d1(self, seq1idx, seq2idx):
        """Given two sequences S and Q, the function d1(S,Q) equals:
        d1(S,Q) = c(SQ)-c(S)+c(QS)-c(Q)
        """
        c1, c2, c12, c21 = self.__get_complexity(seq1idx, seq2idx)
        return c12 - c1 + c21 - c2

    def pwdist_d1_star(self, seq1idx, seq2idx):
        """Given two sequences S and Q, the function d1*(S,Q) equals:
        d1*(S,Q) = [c(SQ)-c(S)+c(QS)-c(Q)]/c(SQ)
        """
        c1, c2, c12, c21 = self.__get_complexity(seq1idx, seq2idx)
        return float(c12 - c1 + c21 - c2) / c12

    def pwdist_d1_star2(self, seq1idx, seq2idx):
        """Given two sequences S and Q, the function d1**(S,Q) equals:
        d1**(S,Q) = \frac{c(SQ)-c(S)+c(QS)-c(Q)}{0.5\cdot [c(SQ)+c(QS)]}
        """

        c1, c2, c12, c21 = self.__get_complexity(seq1idx, seq2idx)
        return float(c12 - c1 + c21 - c2) * 2 / (c12 + c21)

    def set_disttype(self, disttype):
        try:
            pwdist_func = 'self.pwdist_%s' % disttype
            self.pairwise_distance = eval(pwdist_func)
        # Method does not exist.
        except AttributeError:
            msg = 'unknown disttype "%s"' % disttype
            raise ValueError(msg)


def main():
    from .utils import distmatrix
    from .utils.seqrecords import main

    seq_records = main()

    distance = Distance(seq_records)

    l = ['d', 'd_star', 'd1', 'd1_star', 'd1_star2']
    for el in l:
        distance.set_disttype(el)

        matrix = distmatrix.create(seq_records.id_list, distance)
        matrix.display()


if __name__ == '__main__':
    main()
