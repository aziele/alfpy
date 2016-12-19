import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsCommonTest.set_test_data()
        self.script_name = 'create_wordpattern.py'

    def test_arg_word_size_when_no_fasta(self):
        args = ['--word_size', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_word_size_0(self):
        args = ['--fasta', self.filename_pep, '--word_size', '0']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--word_size must be >= 1', out)

    def test_arg_teiresias_when_no_l(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--teiresias']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Teiresias requires --l', out)

    def test_arg_teiresias_when_no_k(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--teiresias', '--l', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Teiresias requires --k', out)

    def test_arg_teiresias_when_k_and_not_l(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--teiresias', '--k', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Teiresias requires --l', out)

    def test_teiresias_when_l_too_small(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--teiresias', '--k', '2', '--l', '1']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--l must be at least 2', out)

    def test_output_word_size_2(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '2aea23ad3e883708dc2f95111f7f04ec')

    def test_output_word_size_2_wordpos(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--word_position']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '040e121be77617191c7d7c847edafc8e')

    def test_output_word_size_1(self):
        args = ['--fasta', self.filename_pep, '--word_size', '1']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '2d4dd98798cb6320975f6919fe43b777')


if __name__ == '__main__':
    unittest.main()
