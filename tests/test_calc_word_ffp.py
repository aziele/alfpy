import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsWordCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsWordCommonTest.set_test_data()
        self.script_name = 'calc_word_ffp.py'

    def test_arg_word_size_when_no_fasta(self):
        args = ['--word_size', '2']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_no_molecule(self):
        args = ['--fasta', self.filename_pep]
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--molecule/-m', out)

    def test_arg_no_word_size(self):
        args = ['--fasta', self.filename_pep, '--molecule', 'protein']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--word_size', out)

    def test_arg_incompatible_args_protein_merge_revcomp(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--molecule', 'protein', '--merge_revcomp']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('Incompatible arguments', out)

    def test_arg_distance_invalid_choice(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--molecule', 'protein', '--distance', 'nonexistent']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('invalid choice', out)

    def test_output_pep_word_size2(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--molecule', 'protein']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '79caa37b67848c52b41a8cb074d810e1')

    def test_output_pep_word_size2_reduce_alphabet(self):
        args = ['--fasta', self.filename_pep, '--word_size', '2',
                '--molecule', 'protein', '--reduce_alphabet']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '2e03fddfa6a10d810c3481fd53ada4a3')

    def test_output_pep_word_pattern2_reduce_alphabet(self):
        args = ['--fasta', self.filename_pep, '--molecule', 'protein',
                '--word_pattern', self.filename_pep_2mer, '--reduce_alphabet']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '2e03fddfa6a10d810c3481fd53ada4a3')

    def test_output_dna_word_size2(self):
        args = ['--fasta', self.filename_dna, '--molecule', 'dna',
                '--word_size', '2']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '69d68abfe5cb8e855f77f9f8fff20178')

    def test_output_dna_word_size2_mergerevcomp(self):
        args = ['--fasta', self.filename_dna, '--molecule', 'dna',
                '--word_size', '2', '--merge_revcomp']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'd3fd336b21aac9922ed7831b8d9f5f83')

    def test_output_dna_word_size2_mergerevcomp_reduce(self):
        args = ['--fasta', self.filename_dna, '--molecule', 'dna',
                '--word_size', '2', '--merge_revcomp', '--reduce_alphabet']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '83fd63884c64c88ee3ff6e4eb2183e8b')

if __name__ == '__main__':
    unittest.main()
