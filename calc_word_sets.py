#! /usr/bin/env python

import argparse
import sys
from alfpy import word_sets_distance
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
    group.add_argument('--word_size', '-s', metavar="N", required=True,
                       help='word size for creating word patterns',
                       type=int)

    group = parser.add_argument_group('OPTIONAL ARGUMENTS')
    distlist = ['dice', 'hamming', 'jaccard']
    group.add_argument('--distance', '-d', choices=distlist,
                       help='choose from: {} [DEFAULT: %(default)s]'.format(
                           ", ".join(distlist)),
                       metavar='', default="dice")

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
        parser.print_usage()
        parser.exit()

    return parser


def validate_args(parser):
    args = parser.parse_args()
    if args.word_size < 1:
        parser.error('Word size must be >= 1.')
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)
    dist = word_sets_distance.Distance(seq_records, args.word_size, args.distance)
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
