#! /usr/bin/env python

# Copyright (c) 2016 Zielezinski A, combio.pl

import argparse
import sys

from alfpy import word_d2
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy.version import __version__


def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate d2 distance between DNA/protein sequences based
        on subsequence (words) occurrences.''',
        add_help=False, prog='calc_word_d2.py'
    )
    group = parser.add_argument_group('REQUIRED ARGUMENTS')
    group.add_argument('--fasta', '-f',
                       help='input FASTA sequence filename', required=True,
                       type=argparse.FileType('r'), metavar="FILE")

    group = parser.add_argument_group('OPTIONAL ARGUMENTS')
    group.add_argument('--min_word_size', '-l',
                       help='minimum word size [default: %(default)s]',
                       type=int, metavar="WORD_SIZE", default=1,
                       )
    group.add_argument('--max_word_size', '-u',
                       help='maximum word size [default: %(default)s]',
                       type=int, metavar="WORD_SIZE", default=3,
                       )
    veclist = ['counts', 'freqs']
    group.add_argument('--vector', '-v', choices=veclist,
                       help='choose from: {} [DEFAULT: %(default)s]'.format(
                            ", ".join(veclist)),
                       metavar='', default="counts")
    group.add_argument('--char_weights', '-W', metavar="FILE",
                       help='''file w/ weights of background sequence characters
                       (nt/aa)''',
                       type=argparse.FileType('r'))

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
    if not args.min_word_size:
        parser.error("min_word_size must be greater than 0")
    elif args.min_word_size >= args.max_word_size:
        parser.error("max_word_size must be greater than min_word_size")
    if args.char_weights:
        try:
            weights = word_vector.read_weightfile(args.char_weights)
            args.char_weights = weights
        except:
            e = 'Invalid format for --char_weights {0}'.format(
                args.char_weights.name)
            parser.error(e)
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)

    patterns = []
    for i in range(args.min_word_size, args.max_word_size + 1):
        p = word_pattern.create(seq_records.seq_list, i)
        patterns.append(p)

    vecs = []
    if args.char_weights is not None:
        weightmodel = word_vector.WeightModel(char_weights=args.char_weights)
        vecklas = {'counts': word_vector.CountsWeight,
                   'freqs': word_vector.FreqsWeight}[args.vector]
        kwargs = {'seq_lengths': seq_records.length_list,
                  'weightmodel': weightmodel}
    else:
        vecklas = {'counts': word_vector.Counts,
                   'freqs': word_vector.Freqs}[args.vector]
        kwargs = {'seq_lengths': seq_records.length_list}
    for p in patterns:
        v = vecklas(patterns=p, **kwargs)
        vecs.append(v)

    dist = word_d2.Distance(vecs)
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
