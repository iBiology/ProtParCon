#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import logging
import unittest

from collections import namedtuple

from Bio import AlignIO

MODEL = namedtuple('model', 'name frequency gamma rates invp type')
PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(PATH, 'tests', 'data', 'utilities')

sys.path.insert(0, PATH)
from imc.utilities import basename, modeling, Tree, trim

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP


class TestBasename(unittest.TestCase):
    result = 'original name'

    def test_fa(self):
        self.assertEqual('a.name', basename('a.name.fa'))

    def test_fan(self):
        self.assertEqual('a.name', basename('a.name.fan'))

    def test_fas(self):
        self.assertEqual('a.name', basename('a.name.fas'))

    def test_fasta(self):
        self.assertEqual('a.name', basename('a.name.fasta'))

    def test_aln(self):
        self.assertEqual('a.name', basename('a.name.aln'))

    def test_clu(self):
        self.assertEqual('a.name', basename('a.name.clu'))

    def test_clustal(self):
        self.assertEqual('a.name', basename('a.name.clustal'))

    def test_phy(self):
        self.assertEqual('a.name', basename('a.name.phy'))

    def test_phylip(self):
        self.assertEqual('a.name', basename('a.name.phylip'))

    def test_phylip_relaxed(self):
        self.assertEqual('a.name', basename('a.name.phylip-relaxed'))

    def test_phy_sequential(self):
        self.assertEqual('a.name', basename('a.name.phylip-sequential'))

    def test_stockholm(self):
        self.assertEqual('a.name', basename('a.name.stockholm'))

    def test_mauve(self):
        self.assertEqual('a.name', basename('a.name.mauve'))

    def test_mau(self):
        self.assertEqual('a.name', basename('a.name.mau'))

    def test_emboss(self):
        self.assertEqual('a.name', basename('a.name.emboss'))

    def test_emb(self):
        self.assertEqual('a.name', basename('a.name.emb'))

    def test_nexus(self):
        self.assertEqual('a.name', basename('a.name.nexus'))

    def test_nex(self):
        self.assertEqual('a.name', basename('a.name.nex'))

    def test_maf(self):
        self.assertEqual('a.name', basename('a.name.maf'))

    def test_xmfa(self):
        self.assertEqual('a.name', basename('a.name.xmfa'))

    def test_newick(self):
        self.assertEqual('a.name', basename('a.name.newick'))

    def test_new(self):
        self.assertEqual('a.name', basename('a.name.new'))

    def test_text(self):
        self.assertEqual('a.name', basename('a.name.text'))

    def test_tsv(self):
        self.assertEqual('a.name', basename('a.name.tsv'))

    def test_name(self):
        self.assertEqual('a.name', basename('a.name'))


class TestModeling(unittest.TestCase):
    def test_model_builtin(self):
        self.assertTupleEqual(MODEL('JTT', 'empirical', 0, 0, 0, 'builtin'),
                              modeling('JTT'))

    def test_model_custom(self):
        self.assertTupleEqual(MODEL((os.path.join(DATA, 'jtt')),
                                    'empirical', 0, 0, 0, 'custom'),
                              modeling(os.path.join(DATA, 'jtt')))

    def test_model_freq(self):
        self.assertTupleEqual(MODEL('JTT', 'empirical', 0, 0, 0, 'builtin'),
                              modeling('JTT+F'))

    def test_model_freq_estimate(self):
        self.assertTupleEqual(MODEL('JTT', 'estimate', 0, 0, 0, 'builtin'),
                              modeling('JTT+FO'))

    def test_model_freq_equal(self):
        self.assertTupleEqual(MODEL('JTT', 'equal', 0, 0, 0, 'builtin'),
                              modeling('JTT+FQ'))

    def test_model_freq_gamma(self):
        self.assertTupleEqual(MODEL('JTT', 'empirical', 4, 0, 0, 'builtin'),
                              modeling('JTT+F+G4'))

    def test_model_freq_estimate_gamma(self):
        self.assertTupleEqual(MODEL('JTT', 'estimate', 4, 0, 0, 'builtin'),
                              modeling('JTT+FO+G4'))

    def test_model_freq_equal_gamma(self):
        self.assertTupleEqual(MODEL('JTT', 'equal', 4, 0, 0, 'builtin'),
                              modeling('JTT+FQ+G4'))

    def test_model_gamma_rate(self):
        self.assertTupleEqual(MODEL('JTT', 'empirical', 4, 8, 0, 'builtin'),
                              modeling('JTT+G4+R8'))

    def test_model_freq_gamma_invp_rate(self):
        self.assertTupleEqual(MODEL('JTT', 'empirical', 4, 4, 'estimate',
                                    'builtin'),
                              modeling('JTT+F+G4+I+R4'))

    def test_model_freq_estimate_gamma_invp_rate(self):
        self.assertTupleEqual(MODEL('JTT', 'estimate', 4, 8, 'estimate',
                                    'builtin'),
                              modeling('JTT+I+FO+G4+R8'))

    def test_model_freq_equal_gamma_invp_rate(self):
        self.assertTupleEqual(MODEL('JTT', 'equal', 4, 4, 'estimate',
                                    'builtin'),
                              modeling('JTT+FQ+G4+I+R4'))

    def test_model_gamma(self):
        self.assertTupleEqual(MODEL('JTT', 'empirical', 4, 0, 0, 'builtin'),
                              modeling('JTT+G4'))

    def test_model_freq_gamma_invp(self):
        self.assertTupleEqual(MODEL('JTT', 'empirical', 4, 0, 'estimate',
                                    'builtin'),
                              modeling('JTT+F+G4+I'))

    def test_model_freq_estimate_gamma_invp(self):
        self.assertTupleEqual(MODEL('JTT', 'estimate', 4, 0, 'estimate',
                                    'builtin'),
                              modeling('JTT+I+FO+G4'))

    def test_model_freq_equal_gamma_invp(self):
        self.assertTupleEqual(MODEL('JTT', 'equal', 4, 0, 'estimate',
                                    'builtin'),
                              modeling('JTT+FQ+G4+I'))


def _read(name):
    with open(name) as f:
        return f.readline().strip()


class TestTree(unittest.TestCase):
    def setUp(self):
        self.top = os.path.join(DATA, 'topology.newick')
        self.top_nodes = os.path.join(DATA, 'topology_internal.newick')
        self.tree = os.path.join(DATA, 'tree.newick')
        self.tree_ic = os.path.join(DATA, 'tree_confidence.newick')
        self.tree_nodes = os.path.join(DATA, 'tree_internal.newick')

        self.s = '(A:0.2,(B:0.3,(C:0.4,''D:0.5):0.2):0.1):0.0;'
        self.s_ic = '(A:0.2,(B:0.3,(C:0.4,D:0.5)7:0.2)6:0.1)5:0.0;'
        self.s_nodes = '(A:0.2,(B:0.3,(C:0.4,D:0.5)E:0.2)F:0.1)root:0.0;'
        
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

    def test_log(self):
        with self.assertLogs('[IMPC]', level='INFO') as cm:
            Tree(self.s + 'abc')
        self.assertRegex(cm.output[0], r'ERROR.*Invalid tree:.*a NEWICK.*.')

    def test_leave(self):
        with self.assertRaises(SystemExit) as cm:
            Tree(self.s + 'abc', leave=True)
        self.assertEqual(cm.exception.code, 1)

    # Test read from file and return a string
    def test_topology_file(self):
        self.assertEqual('(A,(B,(C,D)));',  Tree(self.top).string())

    def test_topology_internal_file(self):
        self.assertEqual('(A,(B,(C,D)));', Tree(self.top_nodes).string())

    def test_tree_file(self):
        self.assertEqual('(A:0.20000,(B:0.30000,(C:0.40000,D:0.50000):0.20000)'
                         ':0.10000):0.00000;', Tree(self.tree).string())

    def test_tree_brlen_file(self):
        self.assertEqual('(A,(B,(C,D)E)F)root;',
                         Tree(self.tree_nodes).string(brlen=False, nodes=True))

    def test_tree_internal_file(self):
        self.assertEqual('(A:0.20000,(B:0.30000,(C:0.40000,D:0.50000)'
                         ':0.20000):0.10000):0.00000;',
                         Tree(self.tree_nodes).string())

    def test_tree_internal_brlen_file(self):
        self.assertEqual('(A,(B,(C,D)));',
                         Tree(self.tree_nodes).string(brlen=False))

    def test_tree_confidence_file(self):
        self.assertEqual('(A,(B,(C,D)7)6)5;',
                         Tree(self.tree_ic).string(ic2name=True, brlen=False))

    # Test read from string and return a string
    def test_topology_string(self):
        self.assertEqual('(A,(B,(C,D)));',
                         Tree('(A,(B,(C,D)));').string())

    def test_topology_internal_string(self):
        self.assertEqual('(A,(B,(C,D)));',
                         Tree('(A,(B,(C,D)E)F)root;').string())

    def test_tree_string(self):
        self.assertEqual('(A:0.20000,(B:0.30000,(C:0.40000,D:0.50000)'
                         ':0.20000):0.10000):0.00000;',
                         Tree(self.s).string())

    def test_tree_brlen_string(self):
        self.assertEqual('(A,(B,(C,D)E)F)root;',
                         Tree(self.s_nodes).string(brlen=False, nodes=True))

    def test_tree_internal_string(self):
        self.assertEqual('(A:0.20000,(B:0.30000,(C:0.40000,D:0.50000):'
                         '0.20000):0.10000):0.00000;',
                         Tree(self.s_nodes).string())

    def test_tree_internal_brlen_string(self):
        self.assertEqual('(A,(B,(C,D)));',
                         Tree(self.s_nodes).string(brlen=False))

    def test_tree_confidence_string(self):
        self.assertEqual('(A,(B,(C,D)7)6)5;',
                         Tree(self.s_ic).string(ic2name=True, brlen=False))

    # Test read from file and write to file
    def test_topology_file_to_file(self):
        name = os.path.join(DATA, 'a.newick')
        self.assertEqual('(A,(B,(C,D)));',
                         _read(Tree(self.top).file(name)))
        self.rm = name

    def test_topology_internal_file_to_file(self):
        name = os.path.join(DATA, 'b.newick')
        self.assertEqual('(A,(B,(C,D)));',
                         _read(Tree(self.top_nodes).file(name)))
        self.rm = name

    def test_tree_file_to_file(self):
        name = os.path.join(DATA, 'c.newick')
        self.assertEqual('(A:0.20000,(B:0.30000,(C:0.40000,D:0.50000):0.20000)'
                         ':0.10000):0.00000;',
                         _read(Tree(self.tree).file(name)))
        self.rm = name

    def test_tree_brlen_file_to_file(self):
        name = os.path.join(DATA, 'd.newick')
        self.assertEqual('(A,(B,(C,D)E)F)root;',
                         _read(Tree(self.tree_nodes).file(name, brlen=False,
                                                          nodes=True)))
        self.rm = name

    def test_tree_internal_file_to_file(self):
        name = os.path.join(DATA, 'e.newick')
        self.assertEqual('(A:0.20000,(B:0.30000,(C:0.40000,D:0.50000)'
                         ':0.20000):0.10000):0.00000;',
                         _read(Tree(self.tree_nodes).file(name)))
        self.rm = name

    def test_tree_internal_brlen_file_to_file(self):
        name = os.path.join(DATA, 'f.newick')
        self.assertEqual('(A,(B,(C,D)));',
                         _read(Tree(self.tree_nodes).file(name, brlen=False)))
        self.rm = name

    def test_tree_confidence_file_to_file(self):
        name = os.path.join(DATA, 'g.newick')
        self.assertEqual('(A,(B,(C,D)7)6)5;',
                         _read(Tree(self.tree_ic).file(name, ic2name=True,
                                                       brlen=False)))
        self.rm = name
        
        
class TestTrim(unittest.TestCase):
    def setUp(self):
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
                
    def test_trim_default(self):
        msa = os.path.join(DATA, 'msa.fa')
        trimmed = {'A': 'ARGPSSSRILAILAVAFIL', 'B': 'ARGPSTSRFLVILAVAFIL',
                   'C': 'ARGPTNSRFLVILAVAFLL', 'D': 'VRVPSTSRFLVILAVAFLL'}
        self.assertDictEqual(trimmed, trim(msa))

    def test_trim_fmt(self):
        msa = os.path.join(DATA, 'msa.phylip')
        trimmed = {'A': 'ARGPSSSRILAILAVAFIL', 'B': 'ARGPSTSRFLVILAVAFIL',
                   'C': 'ARGPTNSRFLVILAVAFLL', 'D': 'VRVPSTSRFLVILAVAFLL'}
        self.assertDictEqual(trimmed, trim(msa, fmt='phylip-relaxed'))

    def test_trim_outfile(self):
        msa = os.path.join(DATA, 'msa.phylip')
        name = os.path.join(DATA, 'trimmed.msa.fa')
        trimmed = {'A': 'ARGPSSSRILAILAVAFIL', 'B': 'ARGPSTSRFLVILAVAFIL',
                   'C': 'ARGPTNSRFLVILAVAFLL', 'D': 'VRVPSTSRFLVILAVAFLL'}
        
        out = trim(msa, fmt='phylip-relaxed', outfile=name)
        self.assertTrue(os.path.isfile(name))
        self.assertDictEqual(trimmed, out)
        rs = AlignIO.read(name, 'fasta')
        out = {r.id: str(r.seq) for r in rs}
        self.assertDictEqual(trimmed, out)
        self.rm = name


if __name__ == '__main__':
    unittest.main()
    pass
