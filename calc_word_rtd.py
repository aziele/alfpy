#! /usr/bin/env python

# Copyright (c) 2016 Zielezinski A, combio.pl

import argparse
import sys

from alfpy import rtd
from alfpy import word_distance
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy.version import __version__

def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate distances between protein/DNA sequences based
        on Return Time Distribution (RTD) of words\' occurrences and their
        relative orders''',
        add_help=False, prog='calc_word_rtd.py'
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
    distlist = word_distance.Distance.get_disttypes()
    group.add_argument('--distance', '-d', choices=distlist,
                       help='choose from: {} [DEFAULT: %(default)s]'.format(
                           ", ".join(distlist)),
                       metavar='', default="google")

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
        if args.word_size < 1:
            parser.error('word size must be >= 1')
    elif args.word_pattern:
        p = word_pattern.read(args.word_pattern)
        if not p.pos_list:
            e = "{0} does not contain info on word positions.\n"
            e += "Please use: create_wordpattern.py with"
            e += " --word_position option."
            parser.error(e.format(args.word_pattern.name))
        else:
            args.word_pattern = p
    else:
        parser.error("Specify either: --word_size or --word_pattern.")
    '''
    try:
        args.alphabet = alphabet.get_alphabet(args.molecule)
    except:
        parser.error("Unknown alphabet {}".format(args.molecule))
    return args
    '''
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)
    if args.word_size:
        p = word_pattern.create(seq_records.seq_list, args.word_size, True)
    else:
        p = args.word_pattern

    vector = rtd.create_vector(seq_records.count, p)
    dist = rtd.Distance(vector, args.distance)

    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
