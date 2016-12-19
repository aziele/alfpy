import hashlib
import os
import subprocess

from alfpy.utils import seqrecords
from alfpy import __version__


ALPHABET_DNA = 'ATGC'
ALPHABET_PEP = 'ACDEFGHIKLMNPRSTQWVY'


def get_test_data(filename):
    filepath = os.path.join(os.path.dirname(__file__), 'data', filename)
    return filepath


def calc_md5(obj):
    return hashlib.md5(str(obj).encode("utf-8")).hexdigest()


def runscript(scriptname, args):
    cmd = ['./' + scriptname]
    for arg in args:
        cmd.append(arg)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True)
    out = "".join(p.communicate())
    return p.returncode, out


class ModulesCommonTest:

    @classmethod
    def set_test_data(cls):
        fh = open(get_test_data('dna.fa'))
        cls.dna_records = seqrecords.read_fasta(fh)
        fh.close()
        fh = open(get_test_data('pep.fa'))
        cls.pep_records = seqrecords.read_fasta(fh)
        fh.close()


class ScriptsCommonTest:
    """Methods testing arguments that are common to all scripts."""

    # the name of the file to read from

    @classmethod
    def set_test_data(cls):
        cls.filename_dna = get_test_data('dna.fa')
        cls.filename_pep = get_test_data('pep.fa')

    def test_arg_version(self):
        cmd = ['--version']
        return_code, out = runscript(self.script_name, cmd)
        self.assertEqual(return_code, 0)
        self.assertIn(__version__, out)

    def test_arg_help(self):
        cmd = ['--help']
        return_code, out = runscript(self.script_name, cmd)
        self.assertEqual(return_code, 0)

    def test_arg_out_when_no_fasta(self):
        cmd = ['--out', 'out.txt']
        return_code, out = runscript(self.script_name, cmd)
        self.assertEqual(return_code, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_outfmt_when_no_fasta(self):
        cmd = ['--outfmt', 'pairwise']
        return_code, out = runscript(self.script_name, cmd)
        self.assertEqual(return_code, 2)
        self.assertIn('--fasta/-f', out)

    def _test_output(self, script_name, args, outfile=True):
        input_filename = args[args.index('--fasta') + 1]
        if outfile:
            args.append('--out')
            output_filename = '{}.out'.format(input_filename)
            args.append(output_filename)
        returncode, result = runscript(script_name, args)
        if outfile:
            fh = open(output_filename)
            result = fh.read()
            fh.close()
            os.remove(output_filename)
        md5 = calc_md5(result)
        return returncode, result, md5


class ScriptsWordCommonTest(ScriptsCommonTest):

    @classmethod
    def set_test_data(cls):
        ScriptsCommonTest.set_test_data()
        cls.filename_char_weights = get_test_data('char_weights.txt')
        cls.filename_char_freqs = get_test_data('char_freqs.txt')
        cls.filename_pep_2mer_wordpos = get_test_data(
            'pep.fa.2mer.wordpos.txt')
        cls.filename_pep_2mer = get_test_data('pep.fa.2mer.txt')
