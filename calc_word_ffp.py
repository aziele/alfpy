#! /usr/bin/env python

# Copyright (c) 2016 Zielezinski A, combio.pl

import argparse
import sys

from alfpy import word_vector
from alfpy import word_distance
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy import word_pattern
from alfpy.utils.data import seqcontent
from alfpy.version import __version__

def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate distance between DNA/protein sequences based
        on feature frequency profiles (FFPs) of words.''',
        add_help=False, prog='calc_word_ffp.py'
    )
    group = parser.add_argument_group('REQUIRED ARGUMENTS')
    group.add_argument('--fasta', '-f',
                       help='input FASTA sequence filename', required=True,
                       type=argparse.FileType('r'), metavar="FILE")
    group.add_argument('--molecule', '-m', choices=['dna', 'rna', 'protein'],
                       help='choose sequence alphabet', required=True)

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
                       metavar='', default="jsd")
    group.add_argument('--reduce_alphabet', '-r', action="store_true",
                       help='''reduce the words' nt/aa alphabet to smaller
                       number of symbols''')
    group.add_argument('--merge_revcomp', '-M', action="store_true",
                       help='''merge together DNA words with their reverse
                       complement words''')

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
        pass
    else:
        parser.error("Specify either: --word_size or --word_pattern.")

    if args.molecule == 'protein' and args.merge_revcomp:
        parser.error("Incompatible arguments: -m protein --merge_revcomp")

    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)
    if args.word_size:
        p = word_pattern.create(seq_records.seq_list, args.word_size)
    else:
        p = word_pattern.read(args.word_pattern)

    if args.reduce_alphabet:
        p = p.reduce_alphabet(seqcontent.get_reduced_alphabet(args.molecule))
    if args.merge_revcomp:
        p.merge_revcomp()

    freqs = word_vector.Freqs(seq_records.length_list, p)

    dist = word_distance.Distance(freqs, args.distance)
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
