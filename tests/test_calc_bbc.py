import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsCommonTest.set_test_data()
        self.script_name = 'calc_bbc.py'

    def test_arg_molecule_when_no_fasta(self):
        args = ['--molecule', 'dna']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_molecule_invalid_choice(self):
        args = ['--fasta', self.filename_dna,
                '--molecule', 'nonexistent_mol']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--molecule/-m', out)

    def test_output_on_dna1(self):
        args = ['--fasta', self.filename_dna, '--m', 'dna']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '6cfc27479ca5fb3d5d2d468544005d8b')

    def test_output_on_dna_k2(self):
        args = ['--fasta', self.filename_dna, '--m', 'dna', '--k', '2']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '1ea7e82d6bb7b8648e0dcca9e089361c')

    def test_output_on_dna_k2_pairwise(self):
        args = ['--fasta', self.filename_dna, '--m', 'dna',
                '--k', '2', '--outfmt', 'pairwise']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '74de6627e68cfb609701c13637ba4090')

    def test_output_on_protein(self):
        args = ['--fasta', self.filename_pep, '--m', 'protein']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '154f2788be2ec349092f22ce359acf80')

    def test_output_on_protein_no_outfile(self):
        args = ['--fasta', self.filename_pep, '--m', 'protein']
        returncode, out, md5 = self._test_output(self.script_name, args, False)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '154f2788be2ec349092f22ce359acf80')


if __name__ == '__main__':
    unittest.main()
