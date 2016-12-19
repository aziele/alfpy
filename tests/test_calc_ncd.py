import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsCommonTest.set_test_data()
        self.script_name = 'calc_ncd.py'

    def test_output_default(self):
        args = ['--fasta', self.filename_pep]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'e5491c3e4197bf1abb92e7f76bdefeaf')

    def test_output_pairwise(self):
        args = ['--fasta', self.filename_pep, '--outfmt', 'pairwise']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'cb69bbabd9a4286a9596f8af3b2b82d5')


if __name__ == '__main__':
    unittest.main()
