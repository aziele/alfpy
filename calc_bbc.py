import argparse
import sys

from alfpy import bbc
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy.utils.data.seqcontent import get_alphabet


def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate distance between DNA/protein sequences based
        on Base-Base Correlation (BBC)''', add_help=False,
    )
    group = parser.add_argument_group('REQUIRED ARGUMENTS')
    group.add_argument('--fasta', '-f',
                       help='input FASTA sequence filename', required=True,
                       type=argparse.FileType('r'), metavar="FILE")
    group.add_argument('--molecule', '-m', choices=['dna', 'rna', 'protein'],
                       help='choose sequence alphabet', required=True)

    group = parser.add_argument_group('OPTIONAL ARGUMENTS')
    group.add_argument('--k', '-k', help='''maximum distance to observe
                        correlation between bases''', type=int, default=10)
    group.add_argument('--out', '-o', help="output filename",
                       metavar="FILE")
    group.add_argument('--outfmt', choices=['phylip', 'pairwise'],
                       default='phylip',
                       help='distances output format [default: %(default)s]')

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
    try:
        args.alphabet = get_alphabet(args.molecule)
    except:
        parser.error("Unknown alphabet {}".format(args.molecule))
    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)
    vector = bbc.create_vectors(seq_records, args.k, alphabet=args.alphabet)
    dist = bbc.Distance(vector)
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
