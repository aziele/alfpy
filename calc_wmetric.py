#! /usr/bin/env python

# Copyright (c) 2016 Zielezinski A, combio.pl

import argparse
import sys

from alfpy import wmetric
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy.utils.data import subsmat
from alfpy.version import __version__


def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate distances between protein sequences based
        on W-metric (Wm).''', add_help=False, prog='calc_wmetric.py'
    )
    group = parser.add_argument_group('REQUIRED ARGUMENTS')
    group.add_argument('--fasta', '-f',
                       help='input FASTA sequence filename', required=True,
                       type=argparse.FileType('r'), metavar="FILE")

    l = subsmat.list_subsmats()
    group = parser.add_argument_group('OPTIONAL ARGUMENTS')
    group.add_argument('--matrix', '-m', choices=l,
                       help='choose from: {} [DEFAULT: %(default)s]'.format(
                           ", ".join(l)), metavar='',
                       default="blosum62")

    group = parser.add_argument_group('OUTPUT ARGUMENTS')
    group.add_argument('--out', '-o', help="output filename",
                       metavar="FILE")
    group.add_argument('--outfmt', choices=['phylip', 'pairwise'],
                       default='phylip',
                       help='distances output format [DEFAULT: %(default)s]')

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
    try:
        args.matrix = subsmat.get(args.matrix)
    except:
        parser.error("Unknown matrix {}".format(args.matrix))
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)
    dist = wmetric.Distance(seq_records, args.matrix)
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
