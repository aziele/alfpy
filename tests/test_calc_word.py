import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsWordCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsWordCommonTest.set_test_data()
        self.script_name = 'calc_word.py'

    def test_arg_word_size_when_no_fasta(self):
        args = ['--word_size', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_word_size_0(self):
        args = ['--fasta', self.filename_pep, '--word_size', '-1']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('error: word size must be >= 1', out)

    def test_arg_no_word_size(self):
        args = ['--fasta', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('error: Specify either: --word_size or', out)

    def test_arg_distance_kls_when_no_freqs(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--distance', 'kld', '--vector', 'counts']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('error: --distance kld requires --vector freqs', out)

    def test_arg_char_weights_when_no_freqs_or_counts(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'freqs_std', '--char_weights',
                self.filename_char_weights]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('error: --char_weights requires a vector of', out)

    def test_arg_freqs_std_when_no_alphabet_size_or_char_freqs(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'freqs_std']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('freqs_std requires either --alphabet_size or ', out)

    def test_arg_freqs_std_when_alphabet_too_small(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'freqs_std', '--alphabet_size', '1']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Alphabet size must be >=2', out)

    def test_arg_char_freqs_when_no_freqs_std(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'counts', '--char_freqs',
                self.filename_char_weights]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--char_freqs requires --vector freqs_std', out)

    def test_arg_alphabet_size_when_no_freqs_std(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'counts', '--alphabet_size', '20']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--alphabet_size requires --vector freqs_std', out)

    def test_arg_char_weights_invalid_format(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'counts',
                '--char_weights', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Invalid format for --char_weights', out)

    def test_arg_char_freqs_invalid_format(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'freqs_std',
                '--char_freqs', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Invalid format for --char_freqs', out)

    def test_arg_distance_invalid_choice(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--distance', 'nonexistent']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('error: argument --distance/-d: invalid choice', out)

    def test_arg_vector_invalid_choice(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'nonexistent']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('error: argument --vector/-v: invalid choice', out)

    def test_output_word_size2_counts_euclid_squared(self):
        # The result of this method is identical to that from decaf+py.
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'counts', '--distance',
                'euclid_squared']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '1b2a44cb139e2ffcb4cc8ebc9eca81a5')

    def test_output_word_size2_freqs_euclid_squared(self):
        # The result of this method is identical to that from decaf+py.
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'freqs', '--distance',
                'euclid_squared']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '043a8af729528379e92618cbeed451ab')

    def test_output_wordsize2_freqs_std_alphabet20_euclid_squared(self):
        # The result of this method is identical to that from decaf+py.
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'freqs_std', '--distance',
                'euclid_squared', '--alphabet_size', '20']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '7cde1e19fd9370e423e8ff3c78dc0562')

    def test_output_wordsize2_freqs_std_char_freqs(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'freqs_std', '--distance',
                'euclid_squared', '--char_freqs', self.filename_char_freqs]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '54d9510aee8c2364bb849a77cab87fce')

    def test_output_word_size2_freqs_euclid_sqaured_pairwise(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--vector', 'freqs', '--distance',
                'euclid_squared', '--outfmt', 'pairwise']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '0f1f15adccf53668a1d2ad776e53bf25')

    def test_output_wordsize2_freqs_char_weights_pairwise(self):
        args = ['--fasta', self.filename_pep,
                '--word_size', '2', '--vector', 'freqs', '--distance',
                'euclid_squared', '--outfmt', 'pairwise',
                '--char_weights', self.filename_char_weights]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'ea1f990dbf28f220496f6a95ff91087b')


if __name__ == '__main__':
    unittest.main()
