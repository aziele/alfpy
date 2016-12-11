#! /usr/bin/env python

# Copyright (c) 2016 Zielezinski A, combio.pl

import argparse
import sys

from alfpy import word_vector
from alfpy import word_distance
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy import word_pattern
from alfpy.version import __version__


def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate compositional distances between DNA/protein
        sequences based on word (of length k) occurrences using a Markov model
        of k-2.''',
        add_help=False, prog='calc_word_comp.py'
    )
    group = parser.add_argument_group('REQUIRED ARGUMENTS')
    group.add_argument('--fasta', '-f',
                       help='input FASTA sequence filename', required=True,
                       type=argparse.FileType('r'), metavar="FILE")

    group = parser.add_argument_group('  Choose between the two options')
    g1 = group.add_mutually_exclusive_group()
    g1.add_argument('--word_size', '-s', metavar="k", type=int,
                    help='''word size (k-mer) for creating word patterns
                        (must be >= 3)'''
                    )
    g1.add_argument('--word_patterns', '-w', nargs=3,
                    help='''3 input word pattern files (k-, [k-1]-,
                        [k-2]-mers)''',
                    type=argparse.FileType('r'), metavar="FILE")

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
    if args.word_size:
        if args.word_size < 3:
            parser.error('Word size must be >= 3')

    elif args.word_patterns:
        l = []
        for i in range(0, 3):
            try:
                p = word_pattern.read(args.word_patterns[i])
                l.append(p)
            except:
                parser.error('Invalid format for word pattern: {0}'.format(
                    args.word_patterns[i].name))

        if len(l) == 3:
            # check if follow rule
            k, k1, k2 = [len(p.pat_list[0]) for p in l]
            if not (k == k1 + 1 == k2 + 2):
                parser.error(
                    '''word pattern lengths do not follow k, k-1, k-2''')

        args.word_patterns = l
    else:
        parser.error("Specify either: --word_size or --word_pattern.")
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)

    if args.word_patterns:
        l = args.word_patterns
    else:
        l = []
        for i in range(args.word_size, args.word_size - 3, -1):
            p = word_pattern.create(seq_records.seq_list, i)
            l.append(p)

    compos = word_vector.Composition(seq_records.length_list, *l)
    dist = word_distance.Distance(compos, 'angle_cos_diss')
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
