#! /usr/bin/env python

import argparse
import sys
from alfpy.utils import seqrecords
from alfpy import word_pattern


def get_parser():
    parser = argparse.ArgumentParser(
        description='''Count subsequences (words) of a given length (size)
        for each sequence in input FASTA-formatted file.''',
        add_help=False,
    )
    group = parser.add_argument_group('REQUIRED ARGUMENTS')
    group.add_argument('--fasta', '-f',
                       help='input FASTA sequence filename', required=True,
                       type=argparse.FileType('r'), metavar="FILE")
    group.add_argument('--word_size', '-w', required=True, type=int,
                       metavar="k", help='word size (>=1)')

    group = parser.add_argument_group('OPTIONAL ARGUMENTS')
    group.add_argument('--word_position', '-p', action="store_true",
                       help='''report word positions in output''')
    group.add_argument('--out', '-o', help="output pattern filename",
                       metavar="FILE")

    t = '  Teiresias options'
    d = '  more info @ https://cm.jefferson.edu/data-tools-downloads/'
    d += 'teiresias-code/\n'
    group = parser.add_argument_group(t, d)
    group.add_argument('--teiresias', '-t', action="store_true",
                       help='''Teiresias program creates word patterns.
                       [by default: disabled]''',
                       )
    group.add_argument('--l', '-l', type=int,
                       help='minimum number of literals and/or brackets')
    group.add_argument('--k', '-k', type=int,
                       help='minimum support that any word can have')

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
    if args.teiresias:
        if args.l is None:
            parser.error("Teiresias requires --l")
        if args.k is None:
            parser.error("Teiresias requires --k")
        if args.word_size < 2:
            parser.error("Teiresias requires --word_size to be >= 2")
        if args.l < 2:
            parser.error("--l must be at least 2")
        if args.l > args.word_size:
            parser.error("--word_size must be >= than --l")
    elif args.word_size < 1:
        parser.error("--word_size must be >= 1")
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    if args.teiresias:
        args.fasta.close()
        p = word_pattern.run_teiresias(args.fasta.name,
                                       w=args.word_size,
                                       l=args.l,
                                       k=args.k,
                                       output_filename=args.out)
    else:
        seq_records = seqrecords.read_fasta(args.fasta)
        args.fasta.close()
        p = word_pattern.create(seq_records.seq_list,
                                args.word_size,
                                args.word_position)

    if args.out:
        oh = open(args.out, 'w')
        oh.write(p.format())
        oh.close()
    else:
        print(p.format())
        # or sys.stdout(p.format()+'\n')


if __name__ == '__main__':
    main()
