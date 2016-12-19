import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsWordCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsWordCommonTest.set_test_data()
        self.script_name = 'calc_word_d2.py'

    def test_arg_when_u_smaller_than_l(self):
        args = ['--fasta', self.filename_pep, '-l', '3', '-u', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('error: max_word_size must be greater than ', out)

    def test_arg_char_weights_invalid_format(self):
        args = ['--fasta', self.filename_pep,
                '-l', '1', '-u', '4',
                '--char_weights', self.filename_pep,
                '--vector', 'freqs']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Invalid format for --char_weights', out)

    def test_arg_word_size_0(self):
        args = ['--fasta', self.filename_pep, '-l', '0']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('min_word_size must be greater than 0', out)

    def test_output_default(self):
        args = ['--fasta', self.filename_pep]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'f651314b77dcd4fe9b3143de28000ca8')

    def test_output_l1_u4(self):
        args = ['--fasta', self.filename_pep, '-l', '1', '-u', '4']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '164ef1a902f74517e6b7cff7798c595f')

    def test_output_l1_u4_freqs(self):
        args = ['--fasta', self.filename_pep, '-l', '1', '-u', '4',
                '--vector', 'freqs']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '8340c1687a0e6ae50c5f6bcc24196247')

    def test_output_l1_u4_char_weights(self):
        args = ['--fasta', self.filename_pep, '-l', '1', '-u', '4',
                '--char_weights', self.filename_char_weights]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '81873a0cb36f7e05698fa664311f38ee')

    def test_script_l1_u4_char_weights_freqs(self):
        args = ['--fasta', self.filename_pep, '-l', '1', '-u', '4',
                '--vector', 'freqs',
                '--char_weights', self.filename_char_weights]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '96c944f9e8e4d2b8ca67bc2620f47d3a')


if __name__ == '__main__':
    unittest.main()
