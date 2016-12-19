import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsWordCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsWordCommonTest.set_test_data()
        self.script_name = 'calc_word_sets.py'

    def test_arg_word_size_when_no_fasta(self):
        args = ['--word_size', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_fasta_when_no_wordsize(self):
        args = ['--fasta', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--word_size', out)

    def test_arg_word_size_too_small(self):
        args = ['--fasta', self.filename_pep, '--word_size', '-1']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Word size must be >= 1.', out)

    def test_arg_distance_invalid_choice(self):
        args = ['--fasta', self.filename_pep, '--word_size', '-1',
                '--distance', 'nonexistent']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('invalid choice', out)

    def test_output_word_size2(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2']
        returncode, out, md5 = self._test_output(self.script_name, args)
        print(out)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'f1b4cf9538d2d2a2a4f1e81ac1b1251d')

    def test_output_word_size2(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--distance', 'jaccard']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '7a744c4665ac06483c5eb36ee03d4fa8')


if __name__ == '__main__':
    unittest.main()
