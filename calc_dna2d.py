#! /usr/bin/env python

# Copyright (c) 2016 Zielezinski A, combio.pl

import argparse
import sys

from alfpy import graph2Ddna
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy.version import __version__


def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate distance between DNA sequences based on
        the two-dimensional (2D) graphical DNA curve''',
        add_help=False, prog='calc_dna2d.py'
    )
    group = parser.add_argument_group('REQUIRED ARGUMENTS')
    group.add_argument('--fasta', '-f',
                       help='input FASTA sequence filename', required=True,
                       type=argparse.FileType('r'), metavar="FILE")

    group = parser.add_argument_group('OPTIONAL ARGUMENTS')
    group.add_argument('--vector', '-v', choices=['2DSV', '2DNV', '2DMV'],
                       help='vector type [default: %(default)s]',
                       default='2DNV')
    group.add_argument('--ndim', '-n', type=int, metavar='N',
                       help='''number of dimensions representing a sequence.
                        (required if --vector 2DMV) [default: %(default)s]''',
                       default=10)

    group = parser.add_argument_group('OUTPUT ARGUMENTS')
    group.add_argument('--out', '-o', help="output filename", metavar="FILE")
    group.add_argument('--outfmt', choices=['phylip', 'pairwise'],
                       default='phylip',
                       help='distances output format [default: %(default)s]')

    group = parser.add_argument_group("OTHER OPTIONS")
    group.add_argument("-h", "--help", action="help",
                       help="show this help message and exit")
    group.add_argument('--version', action='version',
                       version='%(prog)s {}'.format(__version__))

    if len(sys.argv[1:]) == 0:
        # parser.print_help()
        parser.print_usage()
        parser.exit()
    return parser


def validate_args(parser):
    args = parser.parse_args()
    if args.vector == '2DMV' and args.ndim is None:
        parser.error("--vector 2DMV requires the --ndim")
    # TODO: mk as a range
    # stackoverflow.com/questions/18700634/python-argparse-integer-condition-12
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)
    if args.vector == '2DSV':
        vector = graph2Ddna.create_2DSGraphVectors(seq_records)
    elif args.vector == '2DNV':
        vector = graph2Ddna.create_2DNGraphVectors(seq_records)
    else:
        vector = graph2Ddna.create_2DMGraphVectors(seq_records, args.ndim)
    dist = graph2Ddna.Distance(vector)
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
