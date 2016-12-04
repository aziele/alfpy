"""This module reports and handles subsequences of a fixed length (words, 
k-mers) in biological sequences. 

The word patterns are stored in the class `Pattern` that provides information 
about positions of the words in full-length sequences and number of times the 
particular word appears in each sequence.

Word pattern can be retrieved from a set of sequences using either:
    1. the `create` function of this module or 
    2. the `run_teiresias` function that requires the Teiresias 
       (Rigoutsos & Floratos, Bioinformatics 1998) program to be installed 
       https://cm.jefferson.edu/data-tools-downloads/teiresias-code/.

TODO:
    * inspect self.count
    * alfree_format teiresias_format index out of range
"""


class Pattern:
    """Store information on words that are present in sequences. 
    Pattern may include wildcards (dots and brackets).

    Attributes:
        pat_list (list)  : list of words (pattern-strings) 
                             e.g. ['ATGC', 'CGCG', 'GCAT']
        occr_list (list) : number of times every word is present in each seq
                             e.g. [{0:1, 1:2}, {2:2}, {0:3}]
        pos_list (list)  : seq positions of every word in each sequence 
                             e.g. [{0:[234], 1:[200, 204]}, 
                                   {2:[20, 45]}, {0:[3, 8, 34]}]
        count (int)      : number of words

    Note:
        Either `occr_list` or `pos_list` can be an empty list, but not both. 
        `occr_list` is specific for alfree, while `pos_list` is present in 
        the teiresias program.      

    """
    def __init__(self, pat_list, occr_list, pos_list=None):
        """Create a Pattern instance.

        Examples:
            >>> pat = ['ATGC', 'CGCG', 'GCAT']
            >>> occ = [{0:1, 1:2}, {2:2}, {0:3}]
            >>> pos = [{0:[234], 1:[200, 204]}, {2:[20, 45]}, {0:[3, 8, 34]}]
            >>> pattern = Pattern(pat, occ, pos)

        """
        self.pat_list = pat_list     # pattern-strings incl. dots and brackets
        self.pos_list = pos_list if pos_list else []
        self.occr_list = occr_list
        self.count = len(self.pat_list)


    def _alfree_format(self):
        """Return word patterns as a list of 4-element tuples:
            - number of word occurrences in input sequences
            - number of input sequences the word is present
            - word/pattern
            - pairs of integer numbers (seq number: number of words)

        Examples:
            >>> print(p._alfree_format())
            [(3, 2, 'ATGC', '0:1 1:2'), 
             (3, 1, 'GCAT', '0:3'), 
             (2, 1, 'CGCG', '2:2')]

        """
        lines = []
        for i, word in enumerate(self.pat_list):
            occr = self.occr_list[i]
            occr_count = sum(occr.values())
            seqs_count = len(occr.keys())
            l = ['{0}:{1}'.format(n, c) for n,c in occr.items()]
            occr_seqs = " ".join(l)
            lines.append((occr_count, seqs_count, word, occr_seqs))
        lines.sort(reverse=True)
        return lines


    def _teiresias_format(self):
        """Return word patterns as a list of 4-element tuples:
            - number of word occurrences in input sequences
            - number of input sequences the word is present
            - word/pattern
            - pairs of integer numbers (seq number: position in the seq)

        Examples:
            >>> print(p._teiresias_format())
            [(3, 2, 'ATGC', '0 234 1 200 1 204'), 
             (3, 1, 'GCAT', '0 3 0 8 0 34'), 
             (2, 1, 'CGCG', '2 20 2 45')]       

        """
        lines = []
        for i, word in enumerate(self.pat_list):
            offset = self.pos_list[i]
            seqs_count = len(offset.keys())

            l = []
            for seqidx in offset:
                for pos in offset[seqidx]:
                    l.append('{0} {1}'.format(seqidx, pos))
            occr_count = len(l)
            offsets = " ".join(l)                    
            lines.append((occr_count, seqs_count, word, offsets))
        lines.sort(reverse=True)
        return lines


    def format(self, formattype=None):
        """Return Patterns as a string either in Alfree or Teiresias format.

        If formattype is not specified, the format will be auto-detected, 
        depending whether `pos_list` or `occr_list` attributes are empty, 
        Alfree and Teiresias formats, respectively.  

        Examples:
            >>> print(p.format('alfree'))
            3   2   ATGC 0:1 1:2
            3   1   GCAT 0:3
            2   1   CGCG 2:2

            >>> print(p.format('teiresias'))
            3   2   ATGC 0 234 1 200 1 204
            3   1   GCAT 0 3 0 8 0 34
            2   1   CGCG 2 20 2 45

            >>> pat = ['ATGC', 'CGCG', 'GCAT']
            >>> occ = [{0:1, 1:2}, {2:2}, {0:3}]
            >>> pattern = Pattern(pat_list=pat, occr_list=occ, pos_list=[])
            >>> print(pattern.format())
            3   2   ATGC 0:1 1:2
            3   1   GCAT 0:3
            2   1   CGCG 2:2

            >>> pat = ['ATGC', 'CGCG', 'GCAT']
            >>> pos = [{0:[234], 1:[200, 204]}, {2:[20, 45]}, {0:[3, 8, 34]}]
            >>> pattern = Pattern(pat_list=pat, occr_list=[], pos_list=pos)
            >>> print(pattern.format())
            3   2   ATGC 0 234 1 200 1 204
            3   1   GCAT 0 3 0 8 0 34
            2   1   CGCG 2 20 2 45

        """
        func = self._alfree_format
        # Detect format.
        if not formattype:
            if self.pos_list:
                func = self._teiresias_format
        elif formattype=='teiresias':
            func = self._teiresias_format
        lines = []
        for l in func():
            lines.append('{}\t{}\t{} {}'.format(l[0], l[1], l[2], l[3]))
        return "\n".join(lines)


    def reduce_alphabet(self, alphabet_dict):
        """Reduce the words' nt/aa alphabet to smaller number of symbols. 

        For example, four-letter DNA alphabet can be distilled to two-letter 
        purine-pyrimidine encoding and proteins can be represented by 5, 4, 3 
        letters according to their different physical-chemical properties.

        Args:
            alphabet_dict (dict): translation of each sequence symbol to 
               reduced symbol (e.g. {'A': 'R', 'G': 'R', 'T': 'Y', 'C': 'Y'})

        Returns:
            instance of Pattern class

        Examples:
            >>> print(p)
            3   2   ATGC 0:1 1:2
            3   1   GCAT 0:3
            2   1   CGCG 2:2

            >>> alphabet_dict = {'A': 'R', 'C': 'Y', 'T': 'Y', 'G': 'R'}
            >>> p = p.reduce_alphabet(alphabet_dict)
            >>> print(p)
            6   2   RYRY 0:4 1:2
            2   1   YRYR 2:2

        """
        d = {}
        for i,w in enumerate(self.pat_list):
            rw = ''
            for c in w:
                rc = c
                if c in alphabet_dict:
                    rc = alphabet_dict[c]
                rw += rc
            if rw not in d:
                d[rw] = {}
            pos = self.occr_list[i]
            for seqidx in pos:
                if seqidx not in d[rw]:
                    d[rw][seqidx] = 0
                d[rw][seqidx] += pos[seqidx]
        pats = list(d.keys())
        occrs = list(d.values())
        return self.__class__(pat_list=pats, occr_list=occrs, pos_list=[])


    def merge_revcomp(self):
        """Merge together DNA k-mers with their reverse complement words.
        (e.g. A and T in the 1-mer counts are merged together)

        k-mer or its reverse complement are essentially equivalent.

        Returns:
            instance of Pattern

        Examples:
            >>> print(pattern)
            3   2   ATGC 0:1 1:2
            3   1   GCAT 0:3
            2   1   CGCG 2:2
            >>> p = pattern.merge_revcomp()
            >>> print(p)
            6   2   ATGC 0:4 1:2
            2   1   CGCG 2:2

        """
        def revcomp(s):
            d = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R':'Y', 'Y':'R'}
            return "".join(d.get(base, base) for base in s[::-1])

        d = {}
        for i, word in enumerate(self.pat_list):
            revword = revcomp(word)
            word = sorted([word, revword])[0]
            if word not in d:
                d[word] = {}
            pos = self.occr_list[i]
            for seqidx in pos:
                if seqidx not in d[word]:
                    d[word][seqidx] = 0
                d[word][seqidx] += pos[seqidx]
        pats = list(d.keys())
        occrs = list(d.values())
        return self.__class__(pat_list=pats, occr_list=occrs, pos_list=[])


    def __repr__(self):
        return self.format()

    def __str__(self):
        return self.format()




def _create_wordpattern(seq_list, k):
    """Create a word pattern for a given list of sequences and word size.

    Since most of alignment-free distance-calculating algorithms ignore the 
    information about the order of words, this function does not record 
    position of words. Therefore, in a resulting Pattern object, the attribute
    `pos_list` is an empty list.

    Args:
        seq_list (list) : list of sequences
        k (int) : word size

    Returns:
        instance of Pattern

    Examples:
        >>> seqs = ['ATGC', 'CGCG', 'GCAT']
        >>> p = _create_wordpattern(seqs, 1)
        >>> print(p)
        4   3   G 0:1 1:2 2:1
        4   3   C 0:1 1:2 2:1
        2   2   T 0:1 2:1
        2   2   A 0:1 2:1

        >>> p =  _create_wordpattern(seqs, 2)
        >>> print(p)
        3   3   GC 0:1 1:1 2:1
        2   2   AT 0:1 2:1
        2   1   CG 1:2
        1   1   TG 0:1
        1   1   CA 2:1

    """
    d = {}
    for seqidx, seq in enumerate(seq_list):
        for i in range(0, len(seq)-k+1):
            word = seq[i:i+k]
            if word not in d:
                d[word] = {}
            if seqidx not in d[word]:
                d[word][seqidx] = 0
            d[word][seqidx] += 1
    pat_list = list(d.keys())
    occr_list = list(d.values())
    return Pattern(pat_list=pat_list, occr_list=occr_list, pos_list=[])


def _create_wordpattern_positions(seq_list, k):
    """Create a full word pattern for a given list of sequences and word size.

    This function records position of words and the resulting Pattern object
    contains this information.

    Args:
        seq_list (list) : list of sequences
        k (int) : word size

    Returns:
        instance of Pattern

    Examples:
        >>> seqs = ['ATGC', 'CGCG', 'GCAT']
        >>> p = _create_wordpattern_positions(seqs, 1)
        >>> print(p)
        4   3   G 0 2 1 1 1 3 2 0
        4   3   C 0 3 1 0 1 2 2 1
        2   2   T 0 1 2 3
        2   2   A 0 0 2 2

        >>> p =  _create_wordpattern_positions(seqs, 2)
        >>> print(p)
        3   3   GC 0 2 1 1 2 0
        2   2   AT 0 0 2 2
        2   1   CG 1 0 1 2
        1   1   TG 0 1
        1   1   CA 2 1
        
    """
    d = {}
    d1 = {}
    for seqidx, seq in enumerate(seq_list):
        for i in range(0, len(seq)-k+1):
            word = seq[i:i+k]
            if word not in d:
                d[word] = {}
                d1[word] = {}
            if seqidx not in d[word]:
                d[word][seqidx] = []
                d1[word][seqidx] = 0
            d[word][seqidx].append(i)
            d1[word][seqidx]+=1
    pat_list = list(d.keys())
    pos_list = list(d.values())
    occr_list = [d1[w] for w in pat_list]
    return Pattern(pat_list=pat_list, occr_list=occr_list, pos_list=pos_list)


def create(seq_list, word_size=1, wordpos=False):
    """Create a word pattern for a given list of sequences and word size.

    This function can either record or ignore position of words.

    Args:
        seq_list (list) : list of sequences
        k (int) : word size
        wordpos (bool) : record (True) or ignore (False) the word positions

    Returns:
        instance of Pattern

    Examples:
        >>> seqs = ['ATGC', 'CGCG', 'GCAT']
        >>> p = create(seqs, 1)
        >>> print(p)
        4   3   G 0:1 1:2 2:1
        4   3   C 0:1 1:2 2:1
        2   2   T 0:1 2:1
        2   2   A 0:1 2:1

        >>> p = create(seqs, 1, True)
        >>> print(p)        
        4   3   G 0 2 1 1 1 3 2 0
        4   3   C 0 3 1 0 1 2 2 1
        2   2   T 0 1 2 3
        2   2   A 0 0 2 2
        
    """
    if wordpos:
        return _create_wordpattern_positions(seq_list, word_size)
    return _create_wordpattern(seq_list, word_size)


def create_from_fasta(handle, k=1, wordpos=False):
    """Create word patterns (Pattern object) from a FASTA file"""
    seq_records = seqrecords.read_fasta(handle)
    return create(seq_records.seq_list, k=k, wordpos=wordpos)
 

def create_from_bigfasta(filename, k=1):
    """Create word patterns (Pattern object) from a big FASTA file.

    This function does not read full-length sequences into memory, but rather
    reads a file line by line. The function does not record word positions.

    """
    fh = open(filename)
    d = {}
    pat_list = []
    pos_list = []
    len_list = []
    
    extend = ''
    word_idx = -1
    seqnum = -1

    for lineno, line in enumerate(fh):
        line = line.strip()
        if line.startswith('>'):
            word = []
            word_size = 0
            id = line[1:].split()[0]
            seqnum += 1
        else:
            for char in line:
                word.append(char)
                word_size +=1  
                if word_size==k:
                    w = "".join(word)
                    if w not in d:
                        pat_list.append(w)
                        pos_list.append({})
                        len_list.append(len(w))
                        word_idx +=1
                        d[w] = word_idx
                    if seqnum not in pos_list[d[w]]:
                        pos_list[d[w]][seqnum] = 0
                    pos_list[d[w]][seqnum]+=1
                    
                    word.pop(0)
                    word_size = k-1
    fh.close()
    return Pattern(pat_list=pat_list, occr_list=pos_list, pos_list=[])



def run_teiresias(input_filename, w=2, l=2, k=2, output_filename=None):
    """Wrapper function to run Teiresias program.

    Requires Teiresias to be installed.

    Args:
        input_filename (str): input_filename
        l (int): minimum number of literals and/or brackets of the output 
           patterns. Every output pattern will have length at least l. 
        w (int): minimum literal and/or bracket density of the output patterns. 
           In the worst case an output pattern will have length -w and -l 
           literals and/or brackets. The rest w-l will be wildcards. It should
           always be at least equal to -l and followed by a number, e.g. -w6
        k (int): minimum support that any pattern can have. Every reported 
           pattern will have at least -k appearances in the file. It should 
           always be larger than or qual to 2 and followed by a number.

    Returns:
        instance of Pattern

    """
    import subprocess
    output = output_filename
    if not output:
        import uuid
        import os
        output_filename = uuid.uuid4().hex
    cmd = 'teiresias -i{0} -o{1} -l{2} -w{3} -k{4} -c3 -p -s > /dev/null'
    cmd = cmd.format(input_filename, output_filename, l, w, k)
    subprocess.call(cmd, shell=True)
    fh = open(output_filename)
    pattern = read(fh)
    fh.close()
    if not output:
        os.remove(output_filename)
    return pattern


def read(handle):
    """Read word patterns (as a Pattern object) from a file.

    This function autodetects whether the patterns are written either as 
    teiresias or alfree format.

    Returns
        Pattern object.

    """
    fh = handle
    pat_list = [] 
    pos_list = []
    occr_list = []
    for line in fh:
        if not line.startswith('#'):
            l = line.strip().split('\t')[2].split()
            pat = l.pop(0)
            pat_list.append(pat)
            d = {}
            # Alfree
            if ':' in l[0]:
                for el in l:
                    el = el.split(':')
                    seqidx = int(el[0])
                    count = int(el[1])
                    d[seqidx] = count
                occr_list.append(d)
            # Teiresias
            else:   
                for seqidx, pos in zip(l[0::2], l[1::2]):
                    d.setdefault(int(seqidx), []).append(int(pos))
                pos_list.append(d)
                occr_list.append({seqi: len(poss) for seqi, poss in d.items()})
    return Pattern(pat_list=pat_list, pos_list=pos_list, occr_list=occr_list)



def main():
    from .utils import seqrecords
    from .utils.data import seqcontent
    # Manual creation of Pattern instance.
    pat_list = ['ATGC', 'CGCG', 'GCAT']
    occr_list = [{0:1, 1:2}, {2:2}, {0:3}]
    pos_list = [{0:[234], 1:[200, 204]}, {2:[20, 45]}, {0:[3, 8, 34]}]
    pattern = Pattern(pat_list=pat_list, occr_list=occr_list, pos_list=[])

    k = 2

    seq_records = seqrecords.main()


    pattern = create(seq_records.seq_list, k)
    print(pattern)

    p = pattern.reduce_alphabet(seqcontent.get_reduced_alphabet('dna'))
    print(p)

    p = p.merge_revcomp()
    print(p)

    pattern = create(seq_records.seq_list, k, wordpos=True)
    print(pattern)
  

if __name__ == '__main__':
    main()