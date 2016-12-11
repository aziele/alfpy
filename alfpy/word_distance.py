"""Methods for computing distances between 1-D vectors representing sequences.

The following code is inspired by, and was created based on, an excellent
Python code `decaf+py` (http://bioinformatics.org.au/tools/decaf+py/)
originally published in:

    1. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
       doi: 10.1080/10635150701294741.
    2. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
       doi: 10.1080/10635150701294741.

"""

import math
import numpy as np

from .utils import distance


def geometric_mean(values):
    """Calculate geometric mean in numpy array.

    Avoids numerical over- and underflow at the cost of compute time.

    Args:
        values (ndarray)
    Returns:
        float

    """
    inv_len = 1.0 / values.size
    x = values**inv_len
    return x.prod()


class Distance(distance.Distance):
    """Combine vector data with pairwise distance method."""

    def pwdist_euclid_squared(self, seq1idx, seq2idx):
        """Squared Euclidean distance

        References:
            1. Blaisdell BE (1986) Proc Natl Acad Sci U S A 83: 5155-5159.
               doi: 10.1073/pnas.83.14.5155

        """
        return super(Distance, self).pwdist_euclid_squared(seq1idx, seq2idx)

    def pwdist_euclid_norm(self, seq1idx, seq2idx):
        """Euclidean distance

        References:
            1. Vinga & Almeida (2003) Bioinformatics 19(4): 513-523.
               doi: 10.1093/bioinformatics/btg005
            2. http://web.ist.utl.pt/susanavinga/NASC/

        """
        return super(Distance, self).pwdist_euclid_norm(seq1idx, seq2idx)

    def pwdist_euclid_seqlen1(self, seq1idx, seq2idx):
        """A variant of Euclidean distance

        References:
            1. Hohl and Ragan. (2007) Systematic Biology. 56. p. 206-221
               doi: 10.1080/10635150701294741.
            2. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
               doi: 10.1080/10635150701294741.
            3. http://bioinformatics.org.au/tools/decaf+py/

        """
        seqlen = self._vector.seq_lengths
        d = float(self.pwdist_euclid_squared(seq1idx, seq2idx))
        value = d / (float(seqlen[seq1idx] + seqlen[seq2idx]) / 2)
        return value

    def pwdist_euclid_seqlen2(self, seq1idx, seq2idx):
        """A variant of Euclidean distance

        References:
            1. Hohl and Ragan. (2007) Systematic Biology. 56. p. 206-221
               doi: 10.1080/10635150701294741.
            2. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
               doi: 10.1080/10635150701294741.
            3. http://bioinformatics.org.au/tools/decaf+py/

        """
        seqlen = self._vector.seq_lengths
        x = self[seq1idx]
        y = self[seq2idx]
        value = np.sum((x / seqlen[seq1idx]**0.5 -
                        y / seqlen[seq2idx]**0.5)**2)

        return value

    def __angle_cos(self, seq1idx, seq2idx):
        """Cosine of the angle between two vectors in the N-dimensional space
        of composition vectors. The value may vary between -1 and 1.

        Reference:
           1. Berry et al. (1999) SIAM Review, 41(2): p. 335-362
              doi: 10.1137/S0036144598347035

        """
        nom = np.sum(self[seq1idx] * self[seq2idx])
        sum1 = np.sum(self[seq1idx]**2)
        sum2 = np.sum(self[seq2idx]**2)
        value = nom / (math.sqrt(sum1) * math.sqrt(sum2))
        return value

    def pwdist_angle_cos_diss(self, seq1idx, seq2idx):
        """Angled-based composition distance. The distance is normalized
        to the interval (0, 1).

        References:
            1. Hao and Qi (2004) J Bioinform Comput Biol, 2(1): p.1-19
        """
        value = (1 - self.__angle_cos(seq1idx, seq2idx)) / 2
        return value

    def pwdist_angle_cos_evol(self, seq1idx, seq2idx):
        """Angled-based evolutionary distance

        References:
            1. Stuart, Moffett, Leader (2002) Mol Biol Evol 19: 554-562
               doi: 10.1093/oxfordjournals.molbev.a004111
            2. Stuart, Moffett, Baker (2002) Bioinformatics 18: 100-108.
               doi:10.1093/bioinformatics/18.1.100.

        """
        value = -math.log((1 + self.__angle_cos(seq1idx, seq2idx)) / 2)
        return value

    def pwdist_manhattan(self, seq1idx, seq2idx):
        """Manhattan (a.k.a. city block) distance between two vectors."""
        value = np.sum(np.absolute(self[seq1idx] - self[seq2idx]))
        return value

    def pwdist_diff_abs_add(self, seq1idx, seq2idx):
        """
        References:
            1. van Helden J (2004) Bioinformatics 20: 399-406.
               doi: 10.1093/bioinformatics/btg425

        """
        value = np.mean(np.absolute(self[seq1idx] - self[seq2idx]))
        return value

    def pwdist_chebyshev(self, seq1idx, seq2idx):
        """Chebyshev distance between two vectors.

        References:
            1. http://scipy.org/

        """
        u = self[seq1idx]
        v = self[seq2idx]
        return max(abs(u - v))

    def pwdist_braycurtis(self, seq1idx, seq2idx):
        """Bray-Curtis distance between two vectors.

        References:
            1. http://scipy.org/

        """
        u = self[seq1idx]
        v = self[seq2idx]
        return abs(u - v).sum() / abs(u + v).sum()

    def pwdist_diff_abs_mult(self, seq1idx, seq2idx):
        """
        References:
            1. van Helden J (2004) Bioinformatics 20: 399-406.
               doi: 10.1093/bioinformatics/btg425
            2. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
               doi: 10.1080/10635150701294741.
            3. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
               doi: 10.1080/10635150701294741.

        """
        values = np.absolute(self[seq1idx] - self[seq2idx])
        filtered = values[np.nonzero(values)]
        value = geometric_mean(filtered)
        return value

    def pwdist_diff_abs_mult1(self, seq1idx, seq2idx):
        """
        References:
            1. van Helden J (2004) Bioinformatics 20: 399-406.
               doi: 10.1093/bioinformatics/btg425
            2. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
               doi: 10.1080/10635150701294741.
            3. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
               doi: 10.1080/10635150701294741.

        """
        values = np.absolute(self[seq1idx] - self[seq2idx])
        values[values == 0] = 1
        value = geometric_mean(values)
        return value

    def pwdist_diff_abs_mult2(self, seq1idx, seq2idx):
        """
        References:
            1. van Helden J (2004) Bioinformatics 20: 399-406.
               doi: 10.1093/bioinformatics/btg425
            2. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
               doi: 10.1080/10635150701294741.
            3. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
               doi: 10.1080/10635150701294741.

        """
        values = np.absolute(self[seq1idx] - self[seq2idx])
        filtered = values[np.nonzero(values)]
        value = geometric_mean(filtered)
        return value

    def pwdist_kld(self, seq1idx, seq2idx):
        """Kullback-Leibler discrepancy (KL) between two vectors.

        The KL discrepancy between sequences X and Y,
        is computed from their L-tuple (word) frequencies.

        References:
            1. Wu, Hsieh, Li (2001) Biometrics 57: 441-448.
               doi: 10.1111/j.0006-341X.2001.00441.x

        Notes:
            1. KL discrepancy must be computed based on relative
               frequencies (those that sum to 1).
            2. To avoid having an infinite dK L (X, Y) when freqs2 = 0,
               the authors suggest modifying the orifinal formulation
               by adding a unit to both terms of the frequency ratio.

        """
        freqs1 = self[seq1idx] + 1
        freqs2 = self[seq2idx] + 1
        values = freqs1 * np.log2(freqs1 / freqs2)
        value = np.sum(values)
        return value

    def pwdist_jsd(self, seq1idx, seq2idx):
        """Jensen-Shannon Divergence (JSD).

        References:
            1. Sims, Jun, Wu, Kim (2009) Proc Natl Acad Sci USA 106: 2677-2682.
               doi: 10.1073/pnas.0813249106

        """

        x = self[seq1idx]
        y = self[seq2idx]

        # Previous version
        # d1 = x*np.log2(2*x/(x+y))
        # d2 = y*np.log2(2*y/(x+y))
        # d = 0.5*np.sum(d1+d2)
        # return d

        _P = x + 0.0000001
        _Q = y + 0.0000001
        _M = 0.5 * (_P + _Q)

        def entropy(p, q):
            a = np.where(p != 0, p * np.log2(p / q), 0)
            return np.sum(a)

        return 0.5 * (entropy(_P, _M) + entropy(_Q, _M))

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
        return super(Distance, self).pwdist_google(seq1idx, seq2idx)

    def pwdist_lcc(self, seq1idx, seq2idx):
        """Linear correlation coefficient (LCC) between two vectors.
        The Pearson correlation coefficient, which falls between -1 and 1
        is normalized from 0 to 1.

        References:
            1. Petrilli P (1993) Bioinformatics 9:205-209.
               doi: 10.1093/bioinformatics/9.2.205
            2. Vinga S, Almeida J (2003) Bioinformatics 19:513-523.
               doi: 10.1093/bioinformatics/btg005

        """
        x = self[seq1idx]
        y = self[seq2idx]
        n = len(x)
        mx = x.mean()
        my = y.mean()
        xm, ym = x - mx, y - my
        r_num = np.add.reduce(xm * ym)
        r_den = np.sqrt(np.sum(xm * xm, 0) * np.sum(ym * ym, 0))
        r = r_num / r_den

        # Presumably, if abs(r) > 1, then it is only some small artifact
        # of floating point arithmetic.
        r = max(min(r, 1.0), -1.0)
        value = (2 - r - 1) / 2  # Alfree normalization.
        return value

    def pwdist_canberra(self, seq1idx, seq2idx):
        """Compute the Canberra distance between two vectors.

        References:
            1. http://scipy.org/

        Notes:
            When `u[i]` and `v[i]` are 0 for given i, then
            the fraction 0/0 = 0 is used in the calculation.
        """
        u = self[seq1idx]
        v = self[seq2idx]
        d = np.nansum(abs(u - v) / (abs(u) + abs(v)))
        return d

    def pwdist_minkowski(self, seq1idx, seq2idx, p=2):
        """Compute the Minkowski distance between two vectors.

        References:
            1. http://scipy.org/

        """
        u = self[seq1idx]
        v = self[seq2idx]
        if p < 1:
            raise ValueError("p must be at least 1")
        dist = np.linalg.norm(u - v, ord=p)
        return dist


def main():
    from .utils.seqrecords import main as main1
    from .word_vector import main as main2
    from .utils import distmatrix

    seq_records = main1()
    count, freq, freqs_std1, freqs_std2, countw, freqw, compos = main2()

    dist = Distance(compos, 'angle_cos_diss')
    matrix = distmatrix.create(seq_records.id_list, dist)
    matrix.display()

    for dtype in dist.get_disttypes():
        dist = Distance(freq, dtype)
        matrix = distmatrix.create(seq_records.id_list, dist)
        matrix.display()

    print(matrix.data)


if __name__ == '__main__':
    main()
