#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(PATH, 'tests', 'data', 'msa')

sys.path.insert(0, PATH)
from imc.msa import _guess, msa

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP

MUSCLE = settings.MUSCLE
MAFFT = settings.MAFFT
CLUSTAL = settings.CLUSTAL

ALIGNERS = {'muscle': MUSCLE, 'mafft': MAFFT, 'clustal': CLUSTAL}

# It seems use PIPE in Popen will cause a warning:
# ResourceWarning: unclosed file <_io.TextIOWrapper ... >
# No best solution for now, will work on it in future.


class MSATest(unittest.TestCase):
    def setUp(self):
        self.seq = os.path.join(DATA, 'sequence.fa')
        self.rm = ''
        
    def tearDown(self):
        if self.rm and os.path.isfile(self.rm):
            try:
                os.remove(self.rm)
            except OSError:
                logging.warning('MSA testing generated a output file: {} the '
                                'file was not deleted.'.format(self.rm))
                pass
        
    @unittest.skipIf(not all(ALIGNERS.values()), 'No aligners were provided.')
    def test_guess(self):
        aligners = set(ALIGNERS.keys())
        self.assertSetEqual(aligners, 
                            set(_guess(a)[0] for a in ALIGNERS.values()))
                             
    @unittest.skipIf(MUSCLE is None, 'No MUSCLE executable was provided.')
    def test_muscle_default(self):
        aln = msa(MUSCLE, self.seq)
        self.assertTrue(os.path.isfile(aln))
        self.rm = aln

    @unittest.skipIf(MUSCLE is None, 'No MUSCLE executable was provided.')
    def test_muscle_outfile(self):
        out = self.seq.replace('sequence.fa', 'msa.out.muscle.fa')
        aln = msa(MUSCLE, self.seq, outfile=out)
        self.assertTrue(os.path.isfile(aln))
        self.assertEqual(out, aln)
        self.rm = aln
        
    @unittest.skipIf(MAFFT is None, 'No MAFFT executable was provided.')
    def test_mafft_default(self):
        aln = msa(MAFFT, self.seq)
        self.assertTrue(os.path.isfile(aln))
        self.rm = aln

    @unittest.skipIf(MAFFT is None, 'No MAFFT executable was provided.')
    def test_mafft_outfile(self):
        out = self.seq.replace('sequence.fa', 'msa.out.mafft.fa')
        aln = msa(MAFFT, self.seq, outfile=out)
        self.assertTrue(os.path.isfile(aln))
        self.assertEqual(out, aln)
        self.rm = aln
        
    @unittest.skipIf(CLUSTAL is None, 'No CLUSTAL executable was provided.')
    def test_clustal_default(self):
        aln = msa(CLUSTAL, self.seq)
        self.assertTrue(os.path.isfile(aln))
        self.rm = aln

    @unittest.skipIf(CLUSTAL is None, 'No CLUSTAL executable was provided.')
    def test_clustal_outfile(self):
        out = self.seq.replace('sequence.fa', 'msa.out.clustal.fa')
        aln = msa(CLUSTAL, self.seq, outfile=out)
        self.assertTrue(os.path.isfile(aln))
        self.assertEqual(out, aln)
        self.rm = aln


if __name__ == '__main__':
    unittest.main()
