"""Reading and writing FASTA format files"""

from itertools import groupby


class FastaRecord():
    """Object representing a Fasta (aka Pearson) record.

    Attributes:
        seq (str)         : Sequence
        id  (str)         : Sequence identifier
        description (str) : Sequence description
    """

    def __init__(self, seq, id, description=None):
        """Create a FastaRecord.

        Example:
            >>> import Fasta
            >>> record = FastaRecord(seq='MRELEAKAT',
            ...                      id='NP_055309.2',
            ...                      description='TNRC6A')
            >>> print(record)
            >NP_055309.2 TNRC6A
            MRELEAKAT
        """
        self.seq = seq
        self.id = id
        self.description = description

    def __iter__(self):
        """Iterate over the letters in the sequence.

        Example:
            >>> import Fasta
            >>> record = Fasta.read(open('sequence.fasta'))
            >>> for amino_acid in record:
            ...     print(amino_acid)
            M
            R
            E
            L
            E

            This is equivalent to iterating over the sequence directly:
            >>> for amino_acid in record.seq:
            ...     print(amino_acid)
            M
            R
            E
            L
            E
        """
        return iter(self.seq)

    def __contains__(self, char):
        """Implements the 'in' keyword, searches the sequence.

        Example:
            >>> import Fasta
            >>> record = Fasta.read(open('sequence.fasta'))
            >>> print('M' in record)
            True
        """
        return char in self.seq

    def __str__(self):
        """Return the record as a string in the fasta format.

        Example:
            >>> import Fasta
            >>> record = FastaRecord(seq='MRELEAKAT',
            ...                      id='NP_055309.2',
            ...                      description='TNRC6A')
            >>> print(record)
            >NP_055309.2 TNRC6A
            MRELEAKAT
        """
        return self.format(wrap=70)

    def __len__(self):
        """Return the length of the sequence.

        Example:
            >>> import Fasta
            >>> record = Fasta.read(open('sequence.fasta'))
            >>> len(record)
            1240
        """
        return len(self.seq)

    def format(self, wrap=70):
        """Return a formatted Fasta record.

        Example:
            >>> import Fasta
            >>> record = SeqRecord(seq='MRELEAKAT',
                                   id='NP_055309.2',
                                   description='TNRC6A')
            >>> print(record.format())
            >NP_055309.2 TNRC6A
            MRELEAKAT
        """
        header = ">{0}".format(self.id)
        if self.description:
            header += " " + self.description
        header += "\n"
        if wrap:
            wseq = []
            for i in range(0, len(self.seq), wrap):
                wseq.append(self.seq[i:i + wrap])
        return header + "\n".join(wseq)


def parse(handle):
    """
    Generator function to iterate over Fasta records (as FastaRecord objects).

    handle - input file containing fasta sequences.
    """
    faiter = (x[1] for x in groupby(handle, lambda l: l[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()

        id = header.split()[0]
        seq = "".join(s.strip() for s in next(faiter))
        desc = header[len(id):].strip()
        yield FastaRecord(seq, id=id, description=desc)


def read(handle):
    """
    Turns a sequence file into a single FastaRecord.

    EXAMPLE:
    >>> import Fasta
    >>> record = Fasta.read(open('sequence.fasta'))
    >>> print(record.id)
    NP_055309.2
    >>> print(record.seq)
    MRELEAKAT

    If the handle contains no records an exception is raised.
    If the handle contains more than one record, the very first one is read.

    Use the Fasta.parse(handle) function if you want
    to read multiple records from the handle.
    """
    iterator = parse(handle)
    try:
        first = next(iterator)
    except StopIteration:
        first = None
    return first


def to_dict(sequences):
    """
    Turns a Fasta sequence iterator or list into a dictionary.

    - sequences: an iterator that returns FastaRecord objects,
      or simply a list of SeqRecord objects.

    Uses record.id as key.

    If there are duplicate keys, an error is raised.

    EXAMPLE:
    >>> import Fasta
    >>> pdict = Fasta.to_dict(Fasta.parse(open('test.fa')))
    >>> print(sorted(pdict.keys()))
    ['gi|195354411|', 'tr|Q8SY33|']
    >>> print(pdict['tr|Q8SY33|'].description)
    Gawky, isoform A [Drosophila melanogaster]
    >>> len(pdict)
    2

    NOTE:
    This approach is not suitable for very large sets of sequences,
    as all the SeqRecord objects are held in memory.
    """
    d = dict()
    for record in sequences:
        key = record.id
        if key in d:
            raise ValueError("Duplicate key '%s'" % key)
        d[key] = record
    return d


if __name__ == '__main__':
    seqs = ['>seq1 desc1', 'ATGCTGATGATAGATG', 'ATGTAGA',
            '>seq2 desc2', 'ATGCTGCT']
    for seq_record in parse(seqs):
        print(seq_record)
