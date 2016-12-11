"""This module creates and handles distance matrices"""

import itertools
import numpy as np
import sys


def create(id_list, distance):
    """Create a distance matrix (as Matrix object).

    Calculate distance measures between all pairs of sequences.

    Args:
        id_list (list): list of sequence identifiers
        distance (obj): instance of distance.Distance

    Returns:
        Matrix object

    Examples:
        >>> vector
        [[ 3.  6.  4.  1.  3.  4.  3.  0.  1.  1.  6.  4.  5.  0.  3.  4.]
         [ 0.  3.  0.  3.  0.  0.  0.  2.  9.  0.  3.  3.  0.  6.  3.  6.]
         [ 9.  0.  0.  3.  0.  0.  0.  2.  6.  0.  3.  3.  0.  3.  3.  3.]]
        >>> disttype = 'minkowski'
        >>> dist = Distance(vector, disttype)
        >>> id_list = ['seq1', 'seq2', 'seq3']
        >>> matrix = create(id_list, dist)

    """
    size = len(id_list)
    rows = np.zeros([size, size])
    for i, j in itertools.combinations(range(size), 2):
        value = distance.pairwise_distance(i, j)
        rows[i][j] = value
        rows[j][i] = value
    # No need to calculate distances between the same sequences.
    # The distance should be zero.
    '''
    for i in range(size):
        value = distance.pairwise_distance(i, i)
        rows[i][i] = value
    '''
    return Matrix(id_list, rows)


def read_highcharts_matrix(id_list, data):
    """Create a distance matrix from a matrix in Highcharts format.

    Args:
        id_list (list): list of sequence identifiers
        data (list of 4-element tuples)
            e.g. [[0, 1, 0.35, 0.19], [0, 2, 1.0, 0.55], [1, 2, 0.88, 0.48]]

    Returns:
        Matrix object
    """
    size = len(id_list)
    rows = np.zeros([size, size])
    for i, j, norm_value, value in data:
        rows[i][j] = value
        rows[j][i] = value
    return Matrix(id_list, rows)


class Matrix():
    """Distance matrix

    Attributes:
        id_list (list): list of sequence identifiers
        data (ndarray): 2-D array of distance values between pairs of seqs

    """

    def __init__(self, id_list, data):
        """
        Example:
            >>> id_list = ['seq1', 'seq2', 'seq3']
            >>> data
            [[ 0.          0.3531587   0.35509333]
             [ 0.3531587   0.          0.295394  ]
             [ 0.35509333  0.295394    0.        ]]
            >>> matrix = Matrix(id_list, data)

        """
        self.id_list = id_list
        self.data = data

    def normalize(self):
        """Normalize distance values to 0-1 range."""
        self.data /= self.max()

    def __iter__(self):
        """Iterate over a distance matrix."""
        size = self.data.shape[0]
        for i, j in itertools.combinations(range(size), 2):
            yield i, j, self.id_list[i], self.id_list[j], self.data[i][j]

    def writer(self, handle, f, decimal_places):
        """Return a distance matrix as a string in `phylip` or `pairwise`
        formats.

        Args:
            handle : output file / sys.stdout
            f (str): phylip / pairwise
            decimal_places (int): round distance value to decimal places

        """
        if f == 'phylip':
            handle.write("   {0}\n".format(len(self.id_list)))
            for i, line in enumerate(self.data):
                # PHYLIP requires that each sequence identifier
                # is maximum 10 characters long.
                seqid = self.id_list[i][:10]
                l = ['{0:.{1}f}'.format(line[i], decimal_places)
                     for i in range(0, len(line))]
                l.insert(0, '{0: <10}'.format(seqid))
                handle.write(" ".join(l) + '\n')
        elif f == 'pairwise':
            for i, j, seqid1, seqid2, distval in self:
                handle.write("{0}\t{1}\t{2:.{3}f}\n".format(seqid1, seqid2, 
                                                           distval,
                                                           decimal_places))

    def display(self, f="phylip", decimal_places=7):
        """Write a distance matrix to the screen."""
        return self.writer(sys.stdout, f, decimal_places)

    def write_to_file(self, handle, f="phylip", decimal_places=7):
        """Write a distance matrix to a file."""
        return self.writer(handle, f, decimal_places)

    def highcharts(self):
        """Return a distance matrix as a list in the Highcharts format."""
        data = []
        maxval = self.max()
        for i, j, seqid1, seqid2, distval in self:
            data.append([i, j, distval / maxval, distval])
        return data

    def format(self, decimal_places=7):
        lines = ["   {0}\n".format(len(self.id_list))]
        for i, line in enumerate(self.data):
            seqid = self.id_list[i][:10]
            l = ['{0:.{1}f}'.format(line[i], decimal_places)
                     for i in range(0, len(line))]
            l.insert(0, '{0: <10}'.format(seqid))
            lines.append(" ".join(l) + '\n')
        return "".join(lines)

    def min(self):
        """Return minimum distance value in matrix"""
        return np.amin(self.data)

    def max(self):
        """Return maximum distance value in matrix"""
        return np.amax(self.data)

    def __repr__(self):
        return str(self.data)


if __name__ == '__main__':
    id_list = ['seq1', 'seq2', 'seq3']
    l = [[0, 0.3531587, 0.35509333],
         [0.3531587, 0, 0.295394],
         [0.35509333, 0.295394, 0.]
         ]
    data = np.array(l)
    matrix = Matrix(id_list, data)
    print(matrix.format())
    print(matrix.highcharts())
