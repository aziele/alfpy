import argparse
import sys

from alfpy import word_distance
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords


def get_parser():
    parser = argparse.ArgumentParser(
        description='''Calculate distances between DNA/protein sequences based
        on subsequence (words) occurrences.''', add_help=False,
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
    veclist = ['counts', 'freqs', 'freqs_std']
    group.add_argument('--vector', '-v', choices=veclist,
                       help='choose from: {} [DEFAULT: %(default)s]'.format(
                            ", ".join(veclist)),
                       metavar='', default="freqs")
    group.add_argument('--char_weights', '-W', metavar="FILE",
                       help='''file w/ weights of background sequence
                       characters (nt/aa)''',
                       type=argparse.FileType('r'))

    group = parser.add_argument_group('FREQUENCY MODEL ARGUMENTS',
                                      '''  Required for vector \'freqs_std\'.
                                      Specify one of the two options:''')
    group.add_argument('--char_freqs', '-F', metavar="FILE",
                       help='''file w/ frequencies of background sequence
                       characters (nt/aa)''',
                       type=argparse.FileType('r'))
    group.add_argument('--alphabet_size', '-a', metavar="N",
                       help='alphabet size', type=int)

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
    if args.word_size:
        if args.word_size < 1:
            parser.error('word size must be >= 1')
    elif args.word_pattern:
        pass
    else:
        parser.error("Specify either: --word_size or --word_pattern.")

    if args.distance == 'kld' and args.vector != 'freqs':
        parser.error("--distance kld requires --vector freqs.")

    if args.char_weights is not None:
        if args.vector == 'freqs_std':
            e = '''--char_weights requires a vector of either \'freqs\' or
             \'counts\''''
            parser.error(e)
        else:
            try:
                weights = word_vector.read_weightfile(args.char_weights)
                args.char_weights = weights
            except:
                e = 'Invalid format for --char_weights {0}'.format(
                    args.char_weights.name)
                parser.error(e)

    if args.vector == 'freqs_std':
        if args.char_freqs is None and args.alphabet_size is None:
            e = "freqs_std requires either --alphabet_size or --char_freqs"
            parser.error(e)
        elif args.char_freqs is not None:
            try:
                freqs = word_vector.read_freqfile(args.char_freqs)
                args.char_freqs = freqs
            except:
                e = 'Invalid format for --freqs {0}'.format(
                    args.char_freqs.name)
                parser.error(e)
        elif args.alphabet_size < 2:
            parser.error('Alphabet size must be >=2.')
    else:
        if args.char_freqs is not None:
            parser.error("Option --char_freqs requires --vector freqs_std ")
        if args.alphabet_size is not None:
            parser.error("Option --alphabet_size requires --vector freqs_std ")

    return args


def main():
    parser = get_parser()
    args = validate_args(parser)

    seq_records = seqrecords.read_fasta(args.fasta)

    if args.word_size:
        p = word_pattern.create(seq_records.seq_list, args.word_size)
    else:
        p = word_pattern.read(args.word_pattern)

    veccls = {'counts': word_vector.Counts,
              'freqs': word_vector.Freqs}
    vecclsw = {'counts': word_vector.CountsWeight,
               'freqs': word_vector.FreqsWeight
               }

    if args.vector == 'counts' or args.vector == 'freqs':
        if args.char_weights is None:
            vec = veccls[args.vector](seq_records.length_list, p)
        else:
            weightmodel = word_vector.WeightModel(
                char_weights=args.char_weights)
            vec = vecclsw[args.vector](seq_records.length_list, p, weightmodel)
    else:
        if args.alphabet_size:
            freqmodel = word_vector.EqualFreqs(
                alphabet_size=args.alphabet_size)
        else:
            freqmodel = word_vector.EquilibriumFreqs(args.char_freqs)
        vec = word_vector.FreqsStd(seq_records.length_list, p, freqmodel)

    dist = word_distance.Distance(vec, args.distance)
    matrix = distmatrix.create(seq_records.id_list, dist)

    if args.out:
        oh = open(args.out, 'w')
        matrix.write_to_file(oh, args.outfmt)
        oh.close()
    else:
        matrix.display(args.outfmt)


if __name__ == '__main__':
    main()
