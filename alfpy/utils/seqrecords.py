from . import fasta


class SeqRecords:
    """Object representing an ordered collection of sequence records.

    Attributes:
        id_list (list)  : List of sequence record identifiers
        seq_list (list) : List of sequence strings
        count (int)     : Number of sequence records
    """

    def __init__(self, id_list=None, seq_list=None):
        """Create a collection (may be empty) of sequence records.

        Example:
            >>> ids = ['seq1', 'seq2']
            >>> seqs = ['ATGCTG', 'TGCTGATAGTA']
            >>> seq_records = SeqRecords(id_list=ids, seq_list=seqs)
            >>> print seq_records
            SeqRecords (noseqs: 2)
        """
        self.count = 0 if not id_list else len(seq_list)
        self.id_list = id_list if id_list else []
        # Make all sequences uppercased.
        self.seq_list = [s.upper() for s in seq_list] if seq_list else []

    def add(self, id, seq):
        """Add a sequence record to the existing collection.

        Args:
            id (str)  : sequence identifier
            seq (str) : sequence string

        Example:
            >>> seq_record.add("seq3", "TGCTGA")
        """
        self.id_list.append(id)
        self.seq_list.append(seq.upper())
        self.count += 1

    def fasta(self):
        """Return sequence records as a mutli-FASTA string.

        Example:
            >>> ids = ['seq1', 'seq2']
            >>> seqs = ['ATGCTG', 'TGCTGATAGTA']
            >>> seq_records = SeqRecords(id_list=ids, seq_list=seqs)
            >>> print seq_records.fasta()
            >seq1
            ATGCTG
            >seq2
            TGCTGATAGTA
        """
        l = []
        for id, seq in self:
            seq_record = fasta.FastaRecord(seq=seq, id=id)
            l.append(seq_record.format())
        return "\n".join(l)

    @property
    def length_list(self):
        """Return a list of the sequences' length_list"""
        return [len(seq) for seq in self.seq_list]

    def __iter__(self):
        """
        Iterate over sequence records in the collection.

        Example:
            >>> for amino_acid in record:
            ...     print(amino_acid)
            seq1
            ATGCTG
            seq2
            TGCTGATAGTA
        """
        for i in range(self.count):
            id = self.id_list[i]
            seq = self.seq_list[i]
            yield id, seq

    def __len__(self):
        """
        Return the number of sequence records in the collection.

        Example:
            >>> len(seq_records)
            3    
        """
        return len(self.seq_list)

    def __repr__(self):
        return "{0} (noseqs: {1})".format(self.__class__.__name__,
                                          self.count)


def read_fasta(handle):
    """Create a SeqRecords object from Fasta file.

    Args:
        file handle : a file containing Fasta sequences.

    """
    id_list = []
    seq_list = []
    for seq_record in fasta.parse(handle):
        id_list.append(seq_record.id)
        seq_list.append(seq_record.seq)
    return SeqRecords(id_list=id_list, seq_list=seq_list)


def main():
    seq_records = SeqRecords()
    seq_records.add(
        'seq1', 'AACGTACCATTGAACGTACCATTGAACGTACCATTGATGCATGGTAGAT')
    seq_records.add('seq2', 'CTAGGGGACTTATCTAGGGGACTTATCTAGGGGACTTAT')
    seq_records.add('seq3', 'CTAGGGAAAATTCTAGGGAAAATTCTAGGGAAAATT')

    import uuid
    import os
    outfilename = uuid.uuid4().hex
    oh = open(outfilename, 'w')
    oh.write(seq_records.fasta())
    oh.close()

    fh = open(outfilename)
    seq_records = read_fasta(fh)
    fh.close()
    os.remove(outfilename)

    return seq_records


if __name__ == '__main__':
    seq_records = main()
    print(seq_records.fasta())
