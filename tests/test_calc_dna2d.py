import unittest

from . import utils


class ScriptTest(unittest.TestCase, utils.ScriptsCommonTest):

    def __init__(self, *args, **kwargs):
        super(ScriptTest, self).__init__(*args, **kwargs)
        utils.ScriptsCommonTest.set_test_data()
        self.script_name = 'calc_dna2d.py'

    def test_arg_vector_when_no_fasta(self):
        args = ['--vector', '2DSV']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('--fasta/-f', out)

    def test_arg_vector_invalid_choice(self):
        args = ['--fasta', self.filename_dna, '--vector', 'nonexistent']
        returncode, out = utils.runscript(self.script_name, args)
        self.assertEqual(returncode, 2)
        self.assertIn('invalid choice', out)

    def test_output_default(self):
        args = ['--fasta', self.filename_dna]
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '496832ba4841a988a46c81770ee54668')

    def test_output_vector_2DSV(self):
        args = ['--fasta', self.filename_dna, '--vector', '2DSV']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, 'e35a44622d4f0411b26e12e8eedcdb64')

    def test_output_vector_2DMV(self):
        args = ['--fasta', self.filename_dna, '--vector', '2DMV']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '7638015e1c25657cd572071f3b9ae7c4')

    def test_script_output_vector_2DNV_pairwise(self):
        args = ['--fasta', self.filename_dna, '--vector', '2DNV',
                '--outfmt', 'pairwise']
        returncode, out, md5 = self._test_output(self.script_name, args)
        self.assertEqual(returncode, 0)
        self.assertEqual(md5, '2921e374b468b6de81a1c9140681a3b4')


if __name__ == '__main__':
    unittest.main()
