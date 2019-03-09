#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, PATH)
from ProtParCon.mlt import _guess, mlt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP

FASTTREE = settings.FASTTREE
RAXML = settings.RAXML
IQTREE = settings.IQTREE
PHYML = settings.PHYML

MLTS = {'FastTree': FASTTREE, 'IQ-TREE': IQTREE,
        'RAxML': RAXML, 'PhyML': PHYML}

def _newick(filename):
    with open(filename) as f:
        line = f.readline().rstrip()
        return True if line.startswith('(') and line.endswith(';') else False


class TestMLT(unittest.TestCase):
    def setUp(self):
        self.msa = os.path.join(PATH, 'tests', 'data', 'mlt', 'msa.fa')
        self.rm = ''
    
    def tearDown(self):
        if CLEANUP and self.rm and os.path.exists(self.rm):
            try:
                if os.path.isfile(self.rm):
                    os.remove(self.rm)
                else:
                    shutil.rmtree(self.rm)
            except OSError as err:
                logging.warning('Failed to cleanup testing data {} due to:'
                                '\n\t{}'.format(self.rm, err))
      
    @unittest.skipIf(not all(MLTS.values()), 'No ML tree inference program provided.')
    def test__guess(self):
        aligners = set(MLTS.keys())
        self.assertSetEqual(aligners,
                            set(_guess(a)[0] for a in MLTS.values()))

    @unittest.skipIf(FASTTREE is None, 'No FastTree executable was provided.')
    def test_fasttree_dfault(self):
        out = mlt(FASTTREE, self.msa)
        self.assertTrue(os.path.isfile(out))
        self.assertTrue(_newick(out))
        self.rm = out

    @unittest.skipIf(FASTTREE is None, 'No FastTree executable was provided.')
    def test_fasttree_outfile(self):
        out = os.path.join(PATH, 'tests', 'data', 'mlt', 'fastatree.ml.tree.newick')
        outfile = mlt(FASTTREE, self.msa, outfile=out)
        self.assertTrue(os.path.isfile(outfile))
        self.assertEqual(out, outfile)
        self.assertTrue(_newick(out))
        self.rm = out

    @unittest.skipIf(IQTREE is None, 'No IQ-TREE executable was provided.')
    def test_iqtree_dfault(self):
        out = mlt(IQTREE, self.msa)
        self.assertTrue(os.path.isfile(out))
        self.assertTrue(_newick(out))
        self.rm = out

    @unittest.skipIf(IQTREE is None, 'No IQ-TREE executable was provided.')
    def test_iqtree_outfile(self):
        out = os.path.join(PATH, 'tests', 'data', 'mlt', 'iqtree.ml.tree.newick')
        outfile = mlt(IQTREE, self.msa, outfile=out)
        self.assertTrue(os.path.isfile(outfile))
        self.assertEqual(out, outfile)
        self.assertTrue(_newick(out))
        self.rm = out
        
    @unittest.skipIf(RAXML is None, 'No RAxML executable was provided.')
    def test_raxml_dfault(self):
        out = mlt(RAXML, self.msa)
        self.assertTrue(os.path.isfile(out))
        self.assertTrue(_newick(out))
        self.rm = out

    @unittest.skipIf(RAXML is None, 'No RAxML executable was provided.')
    def test_raxml_outfile(self):
        out = os.path.join(PATH, 'tests', 'data', 'mlt', 'raxml.ml.tree.newick')
        outfile = mlt(RAXML, self.msa, outfile=out)
        self.assertTrue(os.path.isfile(outfile))
        self.assertEqual(out, outfile)
        self.assertTrue(_newick(out))
        self.rm = out
        
    @unittest.skipIf(PHYML is None, 'No PhyML executable was provided.')
    def test_phyml_dfault(self):
        out = mlt(PHYML, self.msa)
        self.assertTrue(os.path.isfile(out))
        self.assertTrue(_newick(out))
        self.rm = out

    @unittest.skipIf(PHYML is None, 'No PhyML executable was provided.')
    def test_phyml_outfile(self):
        out = os.path.join(PATH, 'tests', 'data', 'mlt', 'phyml.ml.tree.newick')
        outfile = mlt(PHYML, self.msa, outfile=out)
        self.assertTrue(os.path.isfile(outfile))
        self.assertEqual(out, outfile)
        self.assertTrue(_newick(out))
        self.rm = out
        
  
if __name__ == '__main__':
    unittest.main(warnings='ignore')
