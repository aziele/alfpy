import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsWordCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsWordCommonTest.set_test_data()
        self.script_name = 'calc_word_rtd.py'

    def test_arg_word_size_when_no_fasta(self):
        args = ['--word_size', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_fasta_when_no_word_size(self):
        args = ['--fasta', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Specify either: --word_size or --word_pattern.', out)

    def test_arg_word_pattern_invalid_format(self):
        args = ['--fasta', self.filename_pep,
                '--word_pattern', self.filename_pep_2mer]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('does not contain info on word positions', out)

    def test_arg_distance_invalid_choice(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--distance', 'nonexistent']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('invalid choice', out)

    def test_output_word_size_2(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '1e1a089908495d60275c039272e8e45f')

    def test_output_wordpattern(self):
        args = ['--fasta', self.filename_pep,
                '--word_pattern', self.filename_pep_2mer_wordpos]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '1e1a089908495d60275c039272e8e45f')

    def test_output_word_size_1(self):
        args = ['--fasta', self.filename_pep, '--outfmt', 'pairwise',
                '--word_size', '1']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'b4f581dabfa83b2f1ff4f5d367865711')


if __name__ == '__main__':
    unittest.main()
