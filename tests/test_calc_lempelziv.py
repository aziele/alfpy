import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsCommonTest.set_test_data()
        self.script_name = 'calc_lempelziv.py'

    def test_agr_fasta_when_invalid_distance(self):
        args = ['--fasta', self.filename_dna,
                '--distance', 'nonexistent']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('invalid choice', out)

    def test_agr_distance_when_no_fasta(self):
        args = ['--distance', 'd1']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_output_default(self):
        args = ['--fasta', self.filename_pep]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '89d18a9ac1e573743fa0214c48dde40c')

    def test_output_distance_d(self):
        args = ['--fasta', self.filename_pep, '--distance', 'd']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'c71cb1521d0fc9084eee21c8599785ef')

    def test_output_distance_d_star_pairwise(self):
        args = ['--fasta', self.filename_pep, '--distance', 'd_star',
                '--outfmt', 'pairwise']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '3ed3ca10d198fe4f44ea85134dbcb481')


if __name__ == '__main__':
    unittest.main()
