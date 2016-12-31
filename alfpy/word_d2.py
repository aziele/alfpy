"""This module computes distance between DNA/protein sequences based on
the d2 metric.

References:
    1. Hide, Burke, Davison (1994) J Comput Biol 1:199-215.
       doi: 10.1089/cmb.1994.1.199
    2. Vinga S, Almeida J (2003) Bioinformatics 19:513-523.
       doi: 10.1093/bioinformatics/btg005

"""

import math
import numpy as np


class Distance:

    """Combine a list of vectors with distance function."""

    def __init__(self, vector_list):
        self.vector_list = vector_list
        self.pairwise_distance = self.pwdist_d2

    def pwdist_d2(self, seqidx1, seqidx2):
        d2 = 0
        for vector in self.vector_list:
            d_res = np.sum((vector[seqidx1]-vector[seqidx2])**2)
            d2 += d_res
        return d2

    def pwdist_d2_squareroot(self, seqidx1, seqidx2):
        return math.sqrt(self.pwdist_d2(seqidx1, seqidx2))

    def set_disttype(self, disttype):
        try:
            pwdist_func = getattr(self, 'pwdist_{}'.format(disttype))
            self.pairwise_distance = pwdist_func
        # Method does not exist.
        except AttributeError:
            msg = 'unknown disttype "{}"'.format(disttype)
            raise ValueError(msg)


def main():
    from .utils.seqrecords import main
    from .utils.data import seqcontent
    from .utils import distmatrix
    from . import word_pattern
    from . import word_vector

    seq_records = main()

    patterns = []
    for i in range(1, 5+1):
        p = word_pattern.create(seq_records.seq_list, i)
        patterns.append(p)

    counts = []
    for p in patterns:
        c = word_vector.Counts(seq_records.length_list, p)
        counts.append(c)

    countsweight = []
    weights = seqcontent.get_weights('dna')
    weightmodel = word_vector.WeightModel(weights)
    for p in patterns:
        c = word_vector.CountsWeight(seq_records, p, weightmodel)
        countsweight.append(c)
    dist = Distance(countsweight)
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()


if __name__ == '__main__':
    main()
