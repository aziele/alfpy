'''Create vectors of word occurrences (e.g. counts, frequencies, standarized
frequncies, weighted counts/frequencies) for sequences.

The following code is inspired by, and was created based on, an excellent
Python code `decaf+py` (http://bioinformatics.org.au/tools/decaf+py/)
originally published in:
    1. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
       doi: 10.1080/10635150701294741.
    2. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
       doi: 10.1080/10635150701294741.

'''

import math
import numpy as np


class Counts:
    """Store counts of words (as word_pattern.Pattern object) in given sequence
    records (as seqrecords.SeqRecords object).

    Attributes:
        seq_lengths (list)    : List of sequence lengths
        pat_list (list)       : List of words
        patlen (int)          : Length of words
        data (numpy.ndarray)  : Array of counts (cols) for each sequence (rows)

    """

    def __init__(self, seq_lengths, patterns):
        """Create Counts object.

        Args:
            seq_lengths (list)  : List of sequence lengths
            pattern (obj: word_pattern.Pattern)

        """
        self.seq_lengths = seq_lengths
        self.pat_list = patterns.pat_list
        self.patlen = len(patterns.pat_list[0])
        self.data = self._get_counts_occurrence(len(seq_lengths), patterns)

    @staticmethod
    def _get_counts_occurrence(seq_count, patterns):
        """Create a matrix of word counts for sequences.

        Args:
            seq_count (int)  : number of sequences
            pattern (obj: word_pattern.Pattern)

        """
        data = np.empty((seq_count, patterns.count))
        for seqidx in range(seq_count):
            for patidx in range(patterns.count):
                occr_dict = patterns.occr_list[patidx]
                count = occr_dict.get(seqidx, 0)
                data[seqidx, patidx] = count
        return data

    def __getitem__(self, seqidx):
        """Return word counts for given sequence based on its index

        Args:
            seqidx (int)  : index of a sequence

        """
        return self.data[seqidx]

    def __str__(self):
        """Return the matrix of counts as a string."""
        result = str(self.__class__)
        for data in self.data:
            result += '\n' + str(list(data))
        return result


class Bools(Counts):
    """Store word occurrences in sequences as Booleans (True / False)."""

    def __init__(self, seq_lengths, patterns):
        Counts.__init__(self, seq_lengths, patterns)
        self.data = self.data.astype(bool)  # ndarray of bools.


class Freqs(Counts):
    """Store word frequencies in sequences."""

    def __init__(self, seq_lengths, patterns):
        Counts.__init__(self, seq_lengths, patterns)
        self.data = self.__relative_freqs()  # Calculate freqs from counts.

    def __relative_freqs(self):
        """Calculates word frequencies."""
        for seqidx in range(self.data.shape[0]):
            seqlen = self.seq_lengths[seqidx]
            counts = self.data[seqidx]
            total = seqlen - self.patlen + 1
            self.data[seqidx] = counts / total
        return self.data


class FreqsStd(Freqs):
    """Store standarized word frequencies of sequences.

    References:
        1. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
           doi: 10.1080/10635150701294741.
        2. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
           doi: 10.1080/10635150701294741.
        3. Wu, Burke and Davison, Biometrics (1997), 53: p.1431-1439
           doi: 10.2307/2533509

    """

    def __init__(self, seq_lengths, patterns, freqmodel):
        # Standardization takes care of variance of occurrences
        Freqs.__init__(self, seq_lengths, patterns)
        self.data = self.__standardize_freqs(freqmodel)

    def __overlaps(self, freqmodel):
        result = []
        for pattern in self.pat_list:
            # assumption: all occurrences point to same pattern
            overlap = freqmodel.overlap_capability(pattern)
            result.append(overlap)
        return result

    def __standardize_freqs(self, freqmodel):
        overlaps = self.__overlaps(freqmodel)
        for seqidx in range(self.data.shape[0]):
            freqs = self.data[seqidx]
            seqlen = self.seq_lengths[seqidx]
            for pidx, pat in enumerate(self.pat_list):
                freq = freqs[pidx]
                overlap = overlaps[pidx]
                var = freqmodel.var(pat, seqlen, self.patlen, overlap)
                self.data[seqidx, pidx] = freq / math.sqrt(var)
        return self.data


class WeightModel:
    """Weighting model for words."""

    def __init__(self, char_weights, wtype='content'):
        self.char_weights = char_weights
        try:
            func = 'self.{0}'.format(wtype)
            self.compute = eval(func)
        # method does not exist
        except AttributeError:
            msg = 'unknown weight model "%s"' % wtype
            raise ValueError(msg)

    def content(self, vector, patterns):
        l = []
        for pat in patterns.pat_list:
            value = 1.0
            for symbol in pat.upper():
                value *= self.char_weights.get(symbol, 1.0)
            l.append(value)
        l = np.array(l)
        return vector * l


class CountsWeight(Counts):
    """Store weighted counts according to a given weight model."""

    def __init__(self, seq_lengths, patterns, weightmodel):
        Counts.__init__(self, seq_lengths, patterns)
        self.data = weightmodel.compute(self.data, patterns)


class FreqsWeight(Freqs):
    """Store weighted freqs according to a given weight model."""

    def __init__(self, seq_lengths, patterns, weightmodel):
        Freqs.__init__(self, seq_lengths, patterns)
        self.data = weightmodel.compute(self.data, patterns)


class WordModel:
    """Variance of frequencies of words (k-mers) based on their overlap
    capabilities.

    References:
        1. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
           doi: 10.1080/10635150701294741.
        2. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
           doi: 10.1080/10635150701294741.
        3. Gentleman and Mullin, Biometrics (1989): 45(1), p.35-52.
           doi: 10.2307/2532033

    """

    def probabilities(self, word):
        """Calculate P_L, the probability of a word,
        and all probabilities P_k of its prefixes.

        Needs to be implemented by subclasses
        """
        raise NotImplementedError

    def overlap_capability(self, word):
        """Calculate overlap capability of a given word.

        Overlap capability indicates to what extent the prefix and suffix
        of a word are equal (i.e. if the word beginning is the same as the
        ending).

        In other words, it indicates periodicity in the word, which leads
        to higher probability of co-occurence of words sharing the repeated
        motifs.

        For example, the word AAAA can occur between 0 and 17 times within
        a sequence of length 20. The word ACAC cannot occur more than 9
        tiummes within a sequence of length 20 because it has less overlap
        capability. The word ACGT jas no overlap capability.

        Args:
            word (str): k-mer to get the overlap capability
        Returns:
            A list of binary values. For example, the word ACAC gives overlap
            capability of [0,1,0,1] and the word AAAC has overlap capability
            [0,0,0,1]

        """
        value = []
        length = len(word)
        for i in range(1, length):
            v = 1 if word[0:i] == word[length - i:length] else 0
            value.append(v)
        # i == len(word) - overlap capability is 1 (word overlaps itself)
        value.append(1)
        return value

    def var(self, word, seq_len, word_len=None, overlap_capability=None,
            word_probs=None):
        """Calculate the variance of word frequencies.

        The variance depends on overlap capability and is larger for
        words with a high level of overlap capability.

        Original formula (LaTeX):
        var(X) = np(1-np)+p^{2}(n-L)(n-L+1)+2p\sum_{k=1}^{min(L-1,n-1)}
        (n-k)Q_{L-k}(\frac{1}{4})^{k}

        where X: frequency of occurrence of subsequence, i.e. word
              Q: overlap capability of word
              L: length of word
              M: length of sequence
              n = M-L+1 (maximum possible number of words)
              p = (1/4)^L
        assuming DNA sequence, alphabet size = 4

        In this implementation, the alphabet size is not assumed to be 4.

        For arbitrary symbol frequencies, substitute in original formula:
        P_L for p
        P_k for (1/4)^k

        Args:
            seq_len (int): length of sequence in which word occurs
            word_len (int): length of word (may be pre-computed)
            overlap_capability (list): may be pre-computed
            word_probs (list): may be pre-computed
        Returns:
            variance (float) of word's frequency
        """

        # compute only if necessary
        if word_len is None:
            word_len = len(word)
        if overlap_capability is None:
            overlap_capability = self.overlap_capability(word)
        if word_probs is None:
            word_probs = self.probabilities(word)

        # p = P_L == last element
        p = word_probs[-1]

        max_num = seq_len - word_len + 1
        # should be: min_term = min(word_len-1, max_num-1)
        #       but: sum_term uses range() which stops at min_term-1
        min_term = min(word_len, max_num)

        # should be: overlap_capability[word_len-k] == Q_{L-k}
        #       but: overlap_capability starts at 0, Q at 1
        # should be: word_probs[k] == P_k
        #       but: word_probs starts at 0, P_k at 1

        sum_term = [(max_num - k) * overlap_capability[word_len - k - 1] *
                    word_probs[k - 1] for k in range(1, min_term)]
        np = max_num * p
        n_L = max_num - word_len
        value = np * (1 - np) + pow(p, 2) * (n_L) * \
            (n_L + 1) + 2 * p * sum(sum_term)
        return value


class EqualFreqs(WordModel):
    """Standarized word fequencies with word model that assumes equal
    frequencies for all sequence characters (symbols).

    References:
        1. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
           doi: 10.1080/10635150701294741.
        2. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
           doi: 10.1080/10635150701294741.

    """

    def __init__(self, alphabet_size):
        """Create the EqualFreqs instance.

        Args:
            alphabet_size (int)
        """
        self.alphabet_size = alphabet_size
        # assumption of equal frequencies
        self._avg_symbol_frequency = 1.0 / self.alphabet_size

    def probabilities(self, word):
        result = []
        value = 1.0
        # actual symbols are not needed
        for _ in word:
            value *= self._avg_symbol_frequency
            result.append(value)
        return result


class EquilibriumFreqs(WordModel):
    """Standarized word fequencies with word model that assumes different
    frequencies (in equilibrium) for all sequence characters (symbols).

    References:
        1. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
           doi: 10.1080/10635150701294741.
        2. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
           doi: 10.1080/10635150701294741.
    """

    def __init__(self, equilibrium_frequencies):
        """0.01
        equilibrium_frequencies: dict with key : value pairs denoting
                                 symbol : equilibrium_frequency
        """
        self._equilibrium_frequencies = equilibrium_frequencies

    def probabilities(self, word):
        result = []
        value = 1.0
        for symbol in word.upper():
            value *= self._equilibrium_frequencies.get(symbol, 0.0)
            result.append(value)
        return result


class Composition(Counts):
    """Composition vector (word counts) that subtracts random counts background.

    Composition vector of word_size of `k` is created on the basis of (k-2)th
    order Markov Model from the actual counting results for the (k-1) and (k-2)
    strings.

    References:
        1. Hao and Qi (2004) J Bioinform Comput Biol 2: 1-19.
           doi: 10.1142/S0219720004000442
        2. Hohl and Ragan. (2007) Systematic Biology. 56. 2007. p. 206-221
           doi: 10.1080/10635150701294741.
        3. Hohl, Rigoutsos, Ragan (2006) Evol Bioinform Online 2: 359-375.
           doi: 10.1080/10635150701294741.

    """

    def __init__(self, seq_lengths, patterns, patterns1, patterns2):

        Counts.__init__(self, seq_lengths, patterns)
        self._counts1 = Counts(seq_lengths, patterns1)
        self._counts2 = Counts(seq_lengths, patterns2)

        self.__check_patlen()
        self._lookup1 = {p: i for i, p in enumerate(patterns1.pat_list)}
        self._lookup2 = {p: i for i, p in enumerate(patterns2.pat_list)}

        for seqnum in range(len(seq_lengths)):
            seqlen = self.seq_lengths[seqnum]
            self.__composition(seqnum, seqlen)

    def __check_patlen(self):

        if self.patlen != self._counts1.patlen + 1 or \
           self.patlen != self._counts2.patlen + 2:

            msg = 'pattern lengths do not follow n, n-1, n-2'

            raise ValueError(msg)

    def __markov_chain_freqs(self, seqnum, seqlen, patnum, patlen, pattern):
        len1 = seqlen - patlen + 1
        len2 = seqlen - patlen + 2
        len3 = seqlen - patlen + 3

        f = float(self.data[seqnum][patnum])

        patternL = pattern[:-1]
        patternR = pattern[1:]
        patternLR = pattern[1:-1]
        patnumL = self._lookup1[patternL]
        patnumR = self._lookup1[patternR]
        patnumLR = self._lookup2[patternLR]
        fL = self._counts1[seqnum][patnumL]
        fR = self._counts1[seqnum][patnumR]
        fLR = self._counts2[seqnum][patnumLR]

        if fLR != 0:
            f0 = float(fL) * fR / fLR * (float(len1) * len3 / pow(len2, 2))
        else:
            f0 = 0.0

        return f, f0

    def __composition(self, seqnum, seqlen):
        compos = []
        patlen = self.patlen
        for patnum, pattern in enumerate(self.pat_list):
            try:

                f, f0 = self.__markov_chain_freqs(seqnum, seqlen, patnum,
                                                  patlen, pattern)
                if np.isnan(f0) or f0 == 0:
                    value = 0.0
                else:
                    value = (f - f0) / f0

            except ZeroDivisionError:
                value = 0.0
            compos.append(value)
        self.data[seqnum] = compos


def _read_charval_file(handle):
    """Read sequence character frequencies/weights from file.

    Args:
        handle (file/list): input file

    Returns:
        dict : e.g. {'A': 0.0826, 'Q': 0.0393}

    Sample input format:
        # Multi-line information
        # e.g. Frequencies obtained from SwissProt.
        A   0.0826
        Q   0.0393
        L   0.0965
        S   0.0659
        R   0.0553
        E   0.0674
        K   0.0583
        T   0.0534
        ...
    """
    d = {}
    for line in handle:
        if line.strip() and not line.startswith('#'):
            sl = line.split()
            char = sl[0]
            val = float(sl[1].strip())
            d[char] = val
    return d


read_weightfile = _read_charval_file
read_freqfile = _read_charval_file


def main():
    from . import word_pattern
    from .utils.seqrecords import main
    from .utils.data import seqcontent

    # protein
    seq_records = main()
    p = word_pattern.create(seq_records.seq_list, 2)
    count = Counts(seq_records.length_list, p)
    freq = Freqs(seq_records.length_list, p)
    freqmodel = EqualFreqs(alphabet_size=4)
    freqs_std1 = FreqsStd(seq_records.length_list, p, freqmodel)

    dna_freqs = seqcontent.get_freqs('dna')
    freqmodel = EquilibriumFreqs(dna_freqs)
    freqs_std2 = FreqsStd(seq_records.length_list, p, freqmodel)

    weights = seqcontent.get_weights('dna')
    weightmodel = WeightModel(weights)
    countw = CountsWeight(seq_records.length_list, p, weightmodel)
    freqw = FreqsWeight(seq_records.length_list, p, weightmodel)

    patternsk = word_pattern.create(seq_records.seq_list, 3)
    patternsk1 = word_pattern.create(seq_records.seq_list, 2)
    patternsk2 = word_pattern.create(seq_records.seq_list, 1)
    compos = Composition(seq_records.length_list,
                         patternsk, patternsk1, patternsk2)
    return count, freq, freqs_std1, freqs_std2, countw, freqw, compos


if __name__ == '__main__':
    main()
