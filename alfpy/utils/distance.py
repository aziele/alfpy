"""This module contains a `Distance` class that combines vector
with distance function.

"""

import math
import numpy as np


class Distance(object):
    """Combine sequences-representing 2-D array of vectors
    with a distance function.

    Attributes:
        _vector (ndarray)
        _disttype (str): distance method name
        pairwise_distance (func): distance method

    """

    def __getitem__(self, seqnum):
        return self._vector[seqnum]

    @classmethod
    def get_disttypes(cls):
        """Return a list of available distance function names.

        Returns:
            list of strings
        """
        l = [x[7:] for x, y in cls.__dict__.items() if x.startswith('pwdist')]
        l.sort()
        return l

    def set_disttype(self, disttype):
        try:
            pwdist_func = 'self.pwdist_{}'.format(disttype)
            self.pairwise_distance = eval(pwdist_func)
        # Distance method does not exist.
        except AttributeError:
            msg = 'unknown disttype "{}"'.format(disttype)
            raise ValueError(msg)

    def __init__(self, vector, disttype):
        """Create instance of Distance.

        Args:
            vector (ndarray)
            disttype (str)

        Examples:
        >>> vector
        [[ 3.  6.  4.  1.  3.  4.  3.  0.  1.  1.  6.  4.  5.  0.  3.  4.]
         [ 0.  3.  0.  3.  0.  0.  0.  2.  9.  0.  3.  3.  0.  6.  3.  6.]
         [ 9.  0.  0.  3.  0.  0.  0.  2.  6.  0.  3.  3.  0.  3.  3.  3.]]
        >>> disttype = 'minkowski'
        >>> dist = Distance(vector, disttype)

        """
        self.set_disttype(disttype)
        self._vector = vector
        self._disttype = disttype

    def pwdist_euclid_squared(self, seq1idx, seq2idx):
        """Squared Euclidean distance

        References:
            1. Blaisdell BE (1986) Proc Natl Acad Sci U S A 83: 5155-5159.
               doi: 10.1073/pnas.83.14.5155

        """
        value = np.sum((self[seq1idx] - self[seq2idx])**2)
        return value

    def pwdist_euclid_norm(self, seq1idx, seq2idx):
        """Euclidean distance

        References:
            1. Vinga & Almeida (2003) Bioinformatics 19(4): 513-523.
               doi: 10.1093/bioinformatics/btg005
            2. http://web.ist.utl.pt/susanavinga/NASC/

        """
        value = math.sqrt(self.pwdist_euclid_squared(seq1idx, seq2idx))
        return value

    def pwdist_google(self, seq1idx, seq2idx):
        """Normalized Google Distance (NGD).

        The maximum values for NGD is 1.0, which means two sequences are
        totally not similar to each other, and the minimum values for
        NGD is 0.0. Therefore, the similarity of the two sequences can be
        obtained by NGS = 1 - NGD. Two sequences are treated as two different
        web pages and the each word frequency represents terms found in each
        webpage.

        References:
            1. Lee & Rashid (2008) Information Technology, ITSim 2008.
               doi:10.1109/ITSIM.2008.4631601

        """
        v1 = self[seq1idx]
        v2 = self[seq2idx]

        sumwx = float(np.sum(v1))
        sumwy = float(np.sum(v2))

        summin = float(np.sum(np.minimum(v1, v2)))

        ngd = (max([sumwx, sumwy]) - summin) / \
            ((sumwx + sumwy) - min([sumwx, sumwy]))
        return ngd
