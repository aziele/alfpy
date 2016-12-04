"""Distance methods measuring dissimilarity between sets of words.

These methods are also implemented in numpy and provided in the 
`word_bool_distance` module. However, here are their faster 
implemetations based on python sets.
"""

from .utils import distance


def _getwords(seq, word_size):
    """Return a set of words (of a given size) that are present
    in a given sequence.

    Args:
        seq (str)
        word_size (int): >= 1

    Example:
        >>> seq = 'ATGCGTA'
        >>> print(_getwords(seq, 2))
        set(['GT', 'CG', 'GC', 'AT', 'TG', 'TA'])
        
    """
    s = set([])
    for i in range(0, len(seq)-word_size+1):
        word = seq[i:i+word_size]
        s.add(word)
    return s



class Distance(distance.Distance):
    """Combine vector data with pairwise distance methods that measures
    dissimilarity between sets."""

    def __init__(self, seq_records, word_size, disttype='jaccard'):
        """Create an instance of Distance
     
        Args:
            seq_records (SeqRecords obj)
            word_size (int): >= 1

        """
        self._vector = [_getwords(s, word_size) for s in seq_records.seq_list] 
        self.set_disttype(disttype)

    def pwdist_jaccard(self, seq1idx, seq2idx):
        """Jaccard distance is complementary to the Jaccard coefficient
        and is obtained by subtracting the Jaccard coefficient from 1."""
        s1 = self[seq1idx]
        s2 = self[seq2idx]
        return 1 - len(s1 & s2) / float(len(s1 | s2))

    def pwdist_dice(self, seq1idx, seq2idx):
        """Sorensen-Dice coefficient (Czekanowski's binary index)"""
        s1 = self[seq1idx]
        s2 = self[seq2idx]
        return 1 - (2 * len(s1 & s2) / float(len(s1)+len(s2)))

    def pwdist_hamming(self, seq1idx, seq2idx):
        """Hamming distance measures the number of words which are in either 
        of the sets and not in their intersection.

        """
        s1 = self[seq1idx]
        s2 = self[seq2idx]
        return len(s1.symmetric_difference(s2))



def main():
    from .utils.seqrecords import SeqRecords
    from .utils import distmatrix

    seq_records = SeqRecords()
    seq_records.add('seq1', 'MKSTGWHF')
    seq_records.add('seq2', 'MKSSSSTGWGWG')
    seq_records.add('seq3', 'MKSTLKNGTEQ') 
    dist = Distance(seq_records, 2, 'jaccard')
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()

if __name__ == '__main__':
    main()