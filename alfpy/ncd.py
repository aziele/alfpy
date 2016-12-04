"""Normalized compression distance (NCD)

The NCD is a family of distances parametrized with the compressor Z. 
The better Z is, the closer the NCD approaches the NID, and the better 
the results are.

As described in: 
1. Bennett, Gacs, Ming, Vintanyi, Zurek 
   IEEE Transactions on Information Theory 1998. 44(4):1407-1423
   doi: 10.1109/18.681318

2. Li, Chen, Li, Ma, Vitanyi
   IEEE Transactions on Information Theory 2004. 50(12):3250-3264
   doi: 10.1109/TIT.2004.838101

3. https://en.wikipedia.org/wiki/Normalized_compression_distance

"""
import itertools
import zlib


def complexity(s):
    """Compress string and return the size of the compression."""
    s = s.encode("utf-8") # Python 3 fix.
    compr = zlib.compress(s)
    c = float(len(compr))
    return c


class Distance():

    def __init__(self, seq_records):

        self.seq_records = seq_records
        self._complexity = {}
        self.numseqs = seq_records.count
        # Precomputed complexity for input sequences
        # as well as all pairwise concatenated sequences.
        self._complexity = self.__precompute_complexity()

    def __precompute_complexity(self):
        d = {}
        seqs = self.seq_records.seq_list
        # Complexity for single input sequences.
        for seqidx, seq in enumerate(seqs):
            d[(seqidx,)] = complexity(seq)
        # Complexity for pairwise concatenated sequences.
        for i, j in itertools.combinations(range(self.numseqs), 2):
            seq12 = seqs[i] + seqs[j]
            c12 = complexity(seq12)
            d[(i, j)] = c12
        return d

    def pairwise_distance(self, seq1idx, seq2idx):
        """Compute NCD between two sequences. 

        Formula:
        NCD_Z(x,y) = \frac{Z(xy) - \min \{Z(x),Z(y)\}}{\max \{Z(x),Z(y)\}}.

        where: 
        Z(x) is the binary length of the sequence `x` compressed 
        with compressor Z
        """ 
        zx = self._complexity[(seq1idx,)]
        zy = self._complexity[(seq2idx,)]
        zxy = self._complexity[(seq1idx,seq2idx)]
        return (zxy-min([zx,zy]))/max([zx, zy])




if __name__ == '__main__':
    from .utils import distmatrix
    from .utils.seqrecords import main
    seq_records = main()

    dist = Distance(seq_records)
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display('pairwise')