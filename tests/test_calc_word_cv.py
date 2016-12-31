import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsWordCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsWordCommonTest.set_test_data()
        self.script_name = 'calc_word_cv.py'

    def test_word_size_smaller_than_3(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('error: Word size must be >= 3', out)

    def test_word_pattern_only_one_file(self):
        args = ['--fasta', self.filename_pep, '--word_pattern',
                self.filename_pep_2mer]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('expected 3 argument', out)

    def test_word_pattern_not_follow_rule(self):
        args = ['--fasta', self.filename_pep, '--word_pattern',
                self.filename_pep_2mer, self.filename_pep_2mer,
                self.filename_pep_2mer]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn(' do not follow k, k-1, k-2', out)

    def test_fasta_when_no_word_size_or_pattern(self):
        args = ['--fasta', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Specify either: --word_size or --word_pattern', out)

    def test_output_word_size(self):
        args = ['--fasta', self.filename_pep, '--word_size', '3']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '4fbba77e4f7a64601e7d0cb3b0b6878d')

    def test_output_word_pattern(self):
        args = ['--fasta', self.filename_pep, '--word_patterns',
                self.filename_pep_3mer, self.filename_pep_2mer,
                self.filename_pep_1mer
                ]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '4fbba77e4f7a64601e7d0cb3b0b6878d')


if __name__ == '__main__':
    unittest.main()
