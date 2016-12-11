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

CHAR_FREQS = """A   0.0826
Q   0.0393
L   0.0965
S   0.0659
R   0.0553
E   0.0674
K   0.0583
T   0.0534
N   0.0406
G   0.0708
M   0.0241
W   0.0109
D   0.0546
H   0.0227
F   0.0386
Y   0.0292
C   0.0137
I   0.0594
P   0.0471
V   0.0687
"""

CHAR_WEIGHTS = """A   1.21065375303
C   7.29927007299
E   1.48367952522
D   1.8315018315
G   1.41242937853
F   2.59067357513
I   1.6835016835
H   4.40528634361
K   1.71526586621
M   4.14937759336
L   1.03626943005
N   2.46305418719
Q   2.54452926209
P   2.12314225053
S   1.51745068285
R   1.80831826401
T   1.87265917603
W   9.17431192661
V   1.45560407569
Y   3.42465753425
"""


class TestCalcWord(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestCalcWord, self).__init__(*args, **kwargs)
        self.input_filename = '{}'.format(uuid.uuid4().hex)
        self.output_filename = '{}.out'.format(self.input_filename)
        self.char_weights_filename = 'char_weights.txt'
        self.char_freqs_filename = 'char_freqs.txt'

    def setUp(self):
        oh = open(self.input_filename, 'w')
        oh.write(FASTA)
        oh.close()
        oh = open(self.char_weights_filename, 'w')
        oh.write(CHAR_WEIGHTS)
        oh.close()
        oh = open(self.char_freqs_filename, 'w')
        oh.write(CHAR_FREQS)
        oh.close()

    def tearDown(self):
        os.remove(self.input_filename)
        os.remove(self.char_weights_filename)
        os.remove(self.char_freqs_filename)

    def test_script_arguments1(self):
        cmd = ['calc_word.py', '--version']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 0)
        self.assertIn(__version__, tup[1])

    def test_script_arguments2(self):
        cmd = ['calc_word.py', '--help']
        exit_code = subprocess.call(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        self.assertEqual(exit_code, 0)

    def test_script_arguments3(self):
        cmd = ['calc_word.py', '--word_size', '2']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --fasta/-f is required', tup[1])

    def test_script_arguments4(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '-1']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: word size must be >= 1', tup[1])

    def test_script_arguments5(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: Specify either: --word_size or', tup[1])

    def test_script_arguments6(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--distance', 'kld', '--vector', 'counts']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: --distance kld requires --vector freqs', tup[1])

    def test_script_arguments7(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'freqs_std', '--char_weights',
               self.char_weights_filename]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: --char_weights requires a vector of', tup[1])

    def test_script_arguments8(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'freqs_std']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('freqs_std requires either --alphabet_size or ', tup[1])

    def test_script_arguments9(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'freqs_std',
               '--alphabet_size', '1']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('Alphabet size must be >=2', tup[1])

    def test_script_arguments10(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'counts',
               '--char_freqs', self.char_freqs_filename]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('--char_freqs requires --vector freqs_std', tup[1])

    def test_script_arguments11(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'counts',
               '--alphabet_size', '20']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('--alphabet_size requires --vector freqs_std', tup[1])

    def test_script_arguments12(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'counts',
               '--char_weights', self.input_filename]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('Invalid format for --char_weights', tup[1])

    def test_script_arguments13(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'freqs_std',
               '--char_freqs', self.input_filename]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('Invalid format for --char_freqs', tup[1])

    def test_script_arguments14(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--distance', 'nosuchdistance']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --distance/-d: invalid choice', tup[1])

    def test_script_arguments15(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'nosuchvector']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        tup = p.communicate()
        self.assertEqual(p.returncode, 2)
        self.assertIn('error: argument --vector/-v: invalid choice', tup[1])

    def test_script_output1(self):
        # same as decaf+py
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'counts', '--distance',
               'euclid_squared', '--out', self.output_filename
               ]
        exit_code = subprocess.call(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        self.assertEqual(exit_code, 0)
        fh = open(self.output_filename)
        result = fh.read()
        md5 = hashlib.md5(result).hexdigest()
        fh.close()
        os.remove(self.output_filename)
        self.assertEqual(md5, '1b2a44cb139e2ffcb4cc8ebc9eca81a5')

    def test_script_output2(self):
        # same as decaf+py
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'freqs', '--distance',
               'euclid_squared', '--out', self.output_filename
               ]
        exit_code = subprocess.call(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        self.assertEqual(exit_code, 0)
        fh = open(self.output_filename)
        result = fh.read()
        md5 = hashlib.md5(result).hexdigest()
        fh.close()
        os.remove(self.output_filename)
        self.assertEqual(md5, '043a8af729528379e92618cbeed451ab')

    def test_script_output3(self):
        # same as decaf+py
        # maybe also include it in word_vector or word_distance
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'freqs_std', '--distance',
               'euclid_squared', '--alphabet_size', '20',
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
        self.assertEqual(md5, '7cde1e19fd9370e423e8ff3c78dc0562')

    def test_script_output4(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'freqs_std', '--distance',
               'euclid_squared', '--char_freqs', self.char_freqs_filename,
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
        self.assertEqual(md5, '54d9510aee8c2364bb849a77cab87fce')

    def test_script_output5(self):
        cmd = ['calc_word.py', '--fasta', self.input_filename,
               '--word_size', '2', '--vector', 'freqs', '--distance',
               'euclid_squared', '--outfmt', 'pairwise',
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
        self.assertEqual(md5, '0f1f15adccf53668a1d2ad776e53bf25')


if __name__ == '__main__':
    unittest.main()
