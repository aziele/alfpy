#! /usr/bin/env python

# Copyright (c) 2016 Zielezinski A, combio.pl

import argparse
import sys

from alfpy import word_bool_distance
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords


def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate distances between DNA/protein sequences based
        on boolean 1-D vectors of word counting occurrences.''',
        add_help=False,
    )
    group = parser.add_argument_group('REQUIRED ARGUMENTS')
    group.add_argument('--fasta', '-f',
                       help='input FASTA sequence filename', required=True,
                       type=argparse.FileType('r'), metavar="FILE")

    group = parser.add_argument_group('  Choose between the two options')
    g1 = group.add_mutually_exclusive_group()
    g1.add_argument('--word_size', '-s', metavar="N",
                    help='word size for creating word patterns',
                    type=int)
    g1.add_argument('--word_pattern', '-w',
                    help='input filename w/ pre-computed word patterns',
                    type=argparse.FileType('r'), metavar="FILE")

    group = parser.add_argument_group('OPTIONAL ARGUMENTS')
    distlist = word_bool_distance.Distance.get_disttypes()
    group.add_argument('--distance', '-d', choices=distlist,
                       help='choose from: {} [DEFAULT: %(default)s]'.format(
                           ", ".join(distlist)),
                       metavar='', default="jaccard")

    group = parser.add_argument_group('OUTPUT ARGUMENTS')
    group.add_argument('--out', '-o', help="output filename",
                       metavar="FILE")
    group.add_argument('--outfmt', choices=['phylip', 'pairwise'],
                       default='phylip',
                       help='distances output format [DEFAULT: %(default)s]')

    group = parser.add_argument_group("OTHER OPTIONS")
    group.add_argument("-h", "--help", action="help",
                       help="show this help message and exit")

    if len(sys.argv[1:]) == 0:
        # parser.print_help()
        parser.print_usage()  # for just the usage line
        parser.exit()

    return parser


def validate_args(parser):
    args = parser.parse_args()
    if args.word_size:
        if args.word_size < 1:
            parser.error('Word size must be >= 1.')
    elif args.word_pattern:
        pass
    else:
        parser.error("Specify either: --word_size or --word_pattern.")
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)
    if args.word_size:
        p = word_pattern.create(seq_records.seq_list, args.word_size)
    else:
        p = word_pattern.read(args.word_pattern)

    bools = word_vector.Bools(seq_records.length_list, p)
    dist = word_bool_distance.Distance(bools, 'jaccard')
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
