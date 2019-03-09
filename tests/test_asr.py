#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, PATH)
from ProtParCon.asr import _guess, asr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP

CODEML = settings.CODEML
RAXML = settings.RAXML

ASR = {'codeml': CODEML, 'raxml': RAXML}


class TestASR(unittest.TestCase):
    def setUp(self):
        self.msa = os.path.join(PATH, 'tests', 'data', 'asr', 'msa.fa')
        self.tree = os.path.join(PATH, 'tests', 'data', 'asr', 'tree.newick')
        self.rm = ''
    
    def tearDown(self):
        if CLEANUP and self.rm and os.path.isfile(self.rm):
            try:
                os.remove(self.rm)
            except OSError:
                logging.warning('MSA testing generated a output file: {} the '
                                'file was not deleted.'.format(self.rm))
                pass
      
    @unittest.skipIf(not all(ASR.values()), 'No ASR program provided.')
    def test__guess(self):
        aligners = set(ASR.keys())
        aligner = set(_guess(a)[0] for a in ASR.values() if a)
        self.assertTrue(aligner.issubset(aligners))

    @unittest.skipIf(CODEML is None, 'No CODEML executable was provided.')
    def test_codeml_dfault(self):
        out = asr(CODEML, self.msa, self.tree, 'JTT')
        self.assertTrue(os.path.isfile(out))
        self.rm = out

    @unittest.skipIf(CODEML is None, 'No CODEML executable was provided.')
    def test_codeml_outfile(self):
        out = os.path.join(PATH, 'tests', 'data', 'asr', 'codeml.ancestors.tsv')
        outfile = asr(CODEML, self.msa, self.tree, 'JTT', outfile=out)
        self.assertTrue(os.path.isfile(outfile))
        self.assertEqual(out, outfile)
        self.rm = out

    @unittest.skipIf(RAXML is None, 'No RAxML executable was provided.')
    def test_raxml_default(self):
        out = asr(RAXML, self.msa, self.tree, 'JTT')
        self.assertTrue(os.path.isfile(out))
        self.rm = out

    @unittest.skipIf(RAXML is None, 'No RAxML executable was provided.')
    def test_raxml_outfile(self):
        out = os.path.join(PATH, 'tests', 'data', 'asr', 'raxml.ancestors.tsv')
        outfile = asr(RAXML, self.msa, self.tree, 'JTT', outfile=out)
        self.assertTrue(os.path.isfile(outfile))
        self.assertEqual(out, outfile)
        self.rm = out

  
if __name__ == '__main__':
    unittest.main(warnings='ignore')
