import hashlib
import os
import subprocess
import unittest
import uuid

from alfpy import __version__


FASTA = """>seq1
MEVVIRSANFTDNAKIIIVQLNASVEINCTRPNNYTRKGIRIGPGRAVYAAEEIIGDNTLKQVVTKLRE
>seq2
MVIRSANFTDNAKIIIVQLNASVEINCTRPNNNTRKGIRIGPGRAVYAAEEIIGDIRRAHCNIS
>seq3
MFTDNAKIIIVQLNASVEINCTRPNNNTRKGIHIGPGRAFYATGEIIGDIRQAHCNISGAKW
>seq4
MFTDNAKIIIVQLNASVEINCTRPNNNTR"""


class TestCalcWmetric(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestCalcWmetric, self).__init__(*args, **kwargs)
        self.input_filename = '{}'.format(uuid.uuid4().hex)
        self.output_filename = '{}.out'.format(self.input_filename)

    def setUp(self):
        oh = open(self.input_filename, 'w')
        oh.write(FASTA)
        oh.close()

    def tearDown(self):
        os.remove(self.input_filename)

    def test_script_arguments1(self):
        cmd = ['calc_wmetric.py', '--matrix', 'blosum62']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --fasta/-f is required', tup[1])

    def test_script_arguments2(self):
        cmd = ['calc_wmetric.py', '--version']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 0)
        self.assertIn(__version__, tup[1])

    def test_script_arguments3(self):
        cmd = ['calc_wmetric.py', '--help']
        exit_code = subprocess.call(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        self.assertEqual(exit_code, 0)

    def test_script_arguments4(self):
        cmd = ['calc_wmetric.py', '--matrix', 'nonexistingmatrix']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --matrix/-m: invalid choice', tup[1])

    def test_script_arguments5(self):
        cmd = ['calc_wmetric.py', '--matrix', 'nonexistingmatrix',
               '--fasta', self.input_filename]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --matrix/-m: invalid choice', tup[1])

    def test_script_arguments6(self):
        cmd = ['calc_wmetric.py', '--outfmt', 'pairwise']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --fasta/-f is required', tup[1])

    def test_script_arguments7(self):
        cmd = ['calc_wmetric.py', '--out', self.output_filename]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --fasta/-f is required', tup[1])

    def test_script_arguments8(self):
        cmd = ['calc_wmetric.py', '--outfmt', 'pairwise',
               'matrix', 'blosum62', '--out', self.output_filename]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --fasta/-f is required', tup[1])

    def test_script_output1(self):
        cmd = ['calc_wmetric.py', '--fasta', self.input_filename,
               '--out', self.output_filename
               ]
        exit_code = subprocess.call(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        self.assertEqual(exit_code, 0)
        fh = open(self.output_filename)
        result = fh.read()
        md5 = hashlib.md5(result).hexdigest()
        fh.close()
        os.remove(self.output_filename)
        self.assertEqual(md5, '27ad675a7a2e5c2872a8ab495f2d4494')

    def test_script_output2(self):
        cmd = ['calc_wmetric.py', '--fasta', self.input_filename,
               '--out', self.output_filename, '--outfmt', 'phylip'
               ]
        exit_code = subprocess.call(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        self.assertEqual(exit_code, 0)
        fh = open(self.output_filename)
        result = fh.read()
        md5 = hashlib.md5(result).hexdigest()
        fh.close()
        os.remove(self.output_filename)
        self.assertEqual(md5, '27ad675a7a2e5c2872a8ab495f2d4494')

    def test_script_output3(self):
        cmd = ['calc_wmetric.py', '--fasta', self.input_filename,
               '--out', self.output_filename, '--outfmt', 'pairwise'
               ]
        exit_code = subprocess.call(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        self.assertEqual(exit_code, 0)
        fh = open(self.output_filename)
        result = fh.read()
        md5 = hashlib.md5(result).hexdigest()
        fh.close()
        os.remove(self.output_filename)
        self.assertEqual(md5, '195fb45ed46a80473e1d004b9ce40e94')

    def test_script_output4(self):
        cmd = ['calc_wmetric.py', '--fasta', self.input_filename,
               '--out', self.output_filename, '--outfmt', 'phylip',
               '--matrix', 'pam250'
               ]
        exit_code = subprocess.call(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        self.assertEqual(exit_code, 0)
        fh = open(self.output_filename)
        result = fh.read()
        md5 = hashlib.md5(result).hexdigest()
        fh.close()
        os.remove(self.output_filename)
        self.assertEqual(md5, '217ed91de43b091205add32a673cf8fe')


if __name__ == '__main__':
    unittest.main()
