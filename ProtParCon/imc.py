#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module is designed to provide a common and easy way for identifying
parallel and convergent amino acid replacements in orthologous protein
sequences. In order to make this module for general use, function ``ProtParCon()``
is built on top of other modules to facilitate the identification of
parallel and convergent amino acid replacements using a wide range of
sequence data. Depending on the sequence data, optional parameters and external
programs may be required.
"""

import os
import sys
import shutil
import logging
import tempfile
import argparse
import numpy as np

from io import StringIO
from itertools import combinations
from collections import namedtuple, Counter, defaultdict

from Bio import Phylo, SeqIO, AlignIO
from scipy.stats import poisson

from ProtParCon import msa, asr, aut, sim, utilities

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error

AMINO_ACIDS = 'ARNDCQEGHILKMFPSTWYV'


def _pairing(tree):
    """
    Checking whether two branches are sister branch pair or branch pair sharing
    the same evolutionary path.
    
    :param tree: object, a tree object.
    :return: tuple, a list of branches and a list of branch pairs.
    """
    
    def comparable(tree, b1, b2):
        (p1, t1), (p2, t2) = b1[:2], b2[:2]
        if p1 == p2:
            return False
        else:
            t1a = [a.name for a in tree.get_path(t1)]
            t2a = [a.name for a in tree.get_path(t2)]
            if (t1 in t2a) or (t2 in t1a):
                return False
            else:
                return True
    
    branches = []
    for clade in tree.find_clades():
        for child in clade:
            branches.append([clade.name, child.name])
    
    pairs = [(b1, b2) for (b1, b2) in combinations(branches, 2)
             if comparable(tree, b1, b2)]
    return branches, pairs
    

def _load(tsv):
    """
    Load tree, rates, and data blocks from a tsv file.
    
    :param tsv: str, path to the tsv file stores ancestral states or simulated
        sequences.
    :return: tuple, tree, rates (list) and sequence records (defaultdict).
    """
    
    tree, rates, records = None, [], defaultdict(list)
    with open(tsv) as handle:
        for line in handle:
            blocks = line.strip().split()
            if len(blocks) >= 2:
                if blocks[0] == '#TREE' and blocks[1].endswith(';'):
                    tree = Phylo.read(StringIO(blocks[1]), 'newick')
                elif blocks[0] == '#RATES':
                    rates = [float(i) for i in blocks[1:]]
                else:
                    records[blocks[0]].append(blocks[1])

    size = [len(v) for v in records.values()]
    if size:
        size = size[0]
    else:
        error('Invalid sequence file {}, the file does not have tab '
              'separated lines for sequences.'.format(tsv))
        sys.exit(1)
        
    if tree is None:
        error('Invalid sequence file {}, the file does not have a line stores '
              'labeled tree for internal nodes.'.format(tsv))
        sys.exit(1)
    else:
        names = set([clade.name for clade in tree.find_clades()])
        if records:
            ids = set(records.keys())
            if names != ids:
                error('Sequence name space does not match tree name space.')
                sys.exit(1)
        else:
            error('Invalid sequence file {}, the file does not have '
                  'tab separated lines for sequence blocks.'.format(tsv))
            sys.exit(1)

    return tree, rates, records, size


def _sequencing(sequence, tree, aligner, ancestor, wd, asr_model, verbose):
    """
    Identify the type foe sequence file.
    
    :param sequence: str, path to a sequence data file.
    :param tree: str, path to a NEWICK tree file.
    :return: tuple, sequence, alignment, ancestor, and simulation data file.
    """
    
    if tree:
        utilities.Tree(tree, leave=True)
        with open(sequence) as handle:
            line = handle.readline().strip()
            if line.startswith('>'):
                handle.seek(0)
                records = SeqIO.parse(handle, 'fasta')
                lengths = [len(record.seq) for record in records]
            else:
                error('NEWICK format tree provided, but the sequence file '
                      'was not in in FASTA format.')
                sys.exit(1)

        if len(set(lengths)) == 1:
            alignment = sequence
        else:
            if aligner:
                aler, _ = msa._guess(aligner)
                outfile = ''.join([utilities.basename(sequence),
                                  '.{}.fa'.format(aler)])
                if os.path.isfile(outfile):
                    info('Using pre-existed alignment file')
                    alignment = outfile
                else:
                    alignment = msa.msa(aligner, sequence, verbose=verbose,
                                        outfile=outfile)
            else:
                error('FASTA format sequence file provided, but no alignment '
                      'program was provided.')
                sys.exit(1)
                
        trimmed = ''.join([utilities.basename(alignment), '.trimmed.fa'])
        if os.path.isfile(trimmed):
            info('Using pre-existed trimmed alignment file.')
        else:
            utilities.trim(alignment, outfile=trimmed, verbose=verbose)
        
        if trimmed:
            if ancestor:
                if trimmed.endswith('.trimmed.fa'):
                    name = trimmed.replace('.trimmed.fa', '')
                else:
                    name = trimmed
                
                aser, _ = asr._guess(ancestor)
                outfile = ''.join([utilities.basename(name),
                                   '-{}.tsv'.format(aser)])
                if os.path.isfile(outfile):
                    info('Using pre-existed ancestral states sequence file.')
                    sequence = outfile
                else:
                    sequence = asr.asr(ancestor, trimmed, tree, asr_model,
                                       verbose=verbose,
                                       outfile=outfile)
            else:
                error('No ancestral reconstruction program was provided.')
                sys.exit(1)
        else:
            sys.exit(1)
    
    tree, rate, records, size = _load(sequence)
    return tree, rate, records, size, sequence


def _frequencing(record, site=True):
    if isinstance(record, dict):
        tips = [v for k, v in record.items() if not k.startswith('NODE')]
    else:
        tips = record
    nseq, nsites = float(len(tips)), len(tips[0])
    
    if site:
        freq = []
        for i in range(nsites):
            site = [v[i] for v in tips]
            freq.append([site.count(a) / nseq for a in AMINO_ACIDS])
        freq = np.array(freq)
    else:
        counts = Counter(''.join(tips))
        total = float(sum([v for k, v in counts.items() if k in AMINO_ACIDS]))
        freq = ','.join(['{:.8f}'.format(counts.get(aa, 0) / total)
                        for aa in AMINO_ACIDS])
    return freq


def _prob(tree, rates, record, pos, pair, probs, pi):
    (t1p, t1), (t2p, t2) = pair
    mrca = tree.common_ancestor(t1, t2).name
    ancestor = record[mrca][pos]
    times = np.array([tree.distance(c[0], c[1]) for c in
                      [(mrca, t1p), (t1p, t1), (mrca, t2p), (t2p, t2)]])
    sf = _frequencing(record)
    rate = rates[pos]

    ts = np.around(times * rate * 10000).astype(int)
    anc = np.array([1 if ancestor == AMINO_ACIDS[i] else 0 for i in range(20)])
    u = sf[pos, :] / pi

    u.shape = (1, 20)
    um = u.repeat(20, axis=0)
    pm = probs / 100 * um
    for i in range(20):
        pm[i, i] = 1 - (np.sum(pm[i, :]) - pm[i, i])

    eye = np.eye(20)
    t1pP = np.dot(anc, np.linalg.matrix_power(pm, ts[0]))
    t2pP = np.dot(anc, np.linalg.matrix_power(pm, ts[2]))
    t1p = np.dot(eye, np.linalg.matrix_power(pm, ts[1]))
    t2p = np.dot(eye, np.linalg.matrix_power(pm, ts[3]))
    for i in range(20):
        t1p[i, i], t2p[i, i] = 0, 0
    t1pP.shape, t2pP.shape = (20, 1), (20, 1)

    pc = np.sum(
        np.multiply(np.sum(np.multiply(t2pP, t2p), axis=0, keepdims=True),
                    np.multiply(t1pP, t1p)))
    p = np.sum(np.multiply(np.multiply(t1pP, t2pP), np.multiply(t1p, t2p)))
    c = pc - p
    return p, c


def _pc(tree, rates, records, size, length, probs=None, pi=None):
    branches, pairs = _pairing(tree)
    pars, cons, details = defaultdict(list), defaultdict(list), []

    detail = namedtuple('replacement', 'category position pair r1 r2 dataset')
    for i in range(size):
        record = {k: v[i] for k, v in records.items()}

        for pair in pairs:
            (t1p, t1), (t2p, t2) = pair
            name = '-'.join([t1, t2])
            if size == 1 and i == 0 and rates and probs is not None:
                info('Computing expected probability for {}'.format(name))
            po, pe, co, ce = 0, 0, 0, 0
            for pos in range(length):
                s1p, s1 = record[t1p][pos], record[t1][pos]
                s2p, s2 = record[t2p][pos], record[t2][pos]
                if s1p != s1 and s2p != s2:
                    if s1 == s2:
                        if size == 1 and i == 0:
                            label = 'OBSERVATION'
                        else:
                            label = 'SIMULATION-{:4d}'.format(i + 1)
                        r1, r2 = '{}{}'.format(s1p, s1), '{}{}'.format(s2p, s2)
                        if s1p == s2p:
                            po += 1
                            cat = 'P'
                        else:
                            co += 1
                            cat = 'C'
                        details.append(detail(cat, pos, name, r1, r2, label))
                            
                if size == 1 and i == 0 and rates and probs is not None:
                    p, c = _prob(tree, rates, record, pos, pair, probs, pi)
                    pe += p
                    ce += c
            if rates and probs is not None:
                pars[name].extend([po, pe])
                cons[name].extend([co, ce])
            else:
                pars[name].append(po)
                cons[name].append(co)
    
    return tree, pars, cons, details
    

def _load_matrix(model):
    if model.lower() == 'jtt':
        model = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                             'data', 'jtt')
    else:
        if os.path.isfile(model):
            pass
        else:
            error('Unsupported model for computing expected changes, '
                  'calculation aborted.')
            return None, None
    probs, pi = np.zeros((20, 20)), np.zeros((20, ))
    with open(model) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) == 20 and all([i.replace('.', '').isdigit()
                                          for i in fields]):
                pi = np.array([float(i) for i in fields])
                break

        n = 0
        for line in f:
            fields = line.strip().split()
            if len(fields) == 20 and all([i.replace('.', '').isdigit()
                                          for i in fields]):
                probs[n, :] = [float(field) for field in fields]
                n += 1

    return probs, pi
    
    
def imc(sequence, tree='', aligner='', ancestor='', simulator='', save=False,
        asr_model='JTT', exp_model='', n=100, verbose=False):
    
    """
    Identify molecular parallel and convergent changes.
    
    :param sequence: str, path to the sequence data file. Sequence data file
        here covers a wide range of files and formats:
        
        * sequences: raw protein sequence file, need to be in FASTA format
          and a NEWICK format tree is also required for argument tree.
        * msa: multiple sequence alignment file, need to be in FASTA format
          and a NEWICK format tree is also required for argument tree.
        * ancestors: reconstructed ancestral states file, need to be in tsv
          (tab separated) file, the first line needs to start with #TREE,
          second line needs to be a blank line, and the rest lines in the
          file need to be tab separated sequence name (or ID) and amino
          acid sequences.
        * simulations: simulated sequences, need to be in tsv file, the
          first line needs to start with #TREE, second line needs to be
          a blank line, each dataset need to be separated by a blank line
          and inside each dataset block, each line should consist of tab
          separated sequence name (or ID) and amino acid sequences.
          
    :param tree: str, NEWICK format tree string or tree file. This need to be
        set according to argument sequence. if sequence is raw sequence file or
        MSA file, tree is required for guiding ancestral states reconstruction.
        If sequence is ancestors or simulations, then tree is not necessary.
    :param aligner: str, executable of an alignment program.
    :param ancestor: str, executable of an ancestral states reconstruction
        program.
    :param simulator: str, executable of an sequence simulation program.
    :param toper: str, executable of an topology test program.
    :param asr_model: str, model name or model file for ancestral states
        reconstruction.
    :param aut_model: str, model name or model file for topology test.
    :param exp_model: str, model name or model file for estimate expected
        changes based on simulation or replacement probability manupilation.
    :param n: int, number of datasets (or duplicates) should be simulated.
    :param verbose: bool, invoke verbose or silent process mode,
        default: False, silent mode.
    :return: tuple, a dict object of counts of parallel replacements, a dict
        object of counts of convergent replacements, a list consists of details
        of replacements (namedtuple) and the p-value of AU Test (float or None).
    """

    logger.setLevel(logging.INFO if verbose else logging.ERROR)
    
    if os.path.isfile(sequence):
        sequence = os.path.abspath(sequence)
        wd = os.path.dirname(sequence)
    else:
        error('Invalid sequence {}, sequence is not a file or dose not '
              'exist, exited.'.format(sequence))
        sys.exit(1)

    tree, rates, records, size, sequence = _sequencing(sequence, tree, aligner,
                                                       ancestor, wd, asr_model,
                                                       verbose)

    pars, cons, details, aup = None, None, None, None
    h1 = ['Category', 'BranchPair']
    h2 = ['Category', 'Position', 'BranchPair', 'R1', 'R2', 'Dataset']
    
    probs, pi = None, None
    if size == 1:
        h1.append('OBS')
        if exp_model:
            if simulator:
                h1.append('EXP')
                h1.extend(['SIM-{}'.format(i + 1) for i in range(n)])
                
            else:
                probs, pi = _load_matrix(exp_model)
                if probs is not None:
                    h1.append('EXP')
    else:
        h1.append('EXP')
        h1.extend(['SIM-{}'.format(i + 1) for i in range(size)])

    tips = [v[0] for k, v in records.items() if not k.startswith('NODE')]
    length = len(tips[0])
    
    tree, pars, cons, details = _pc(tree, rates, records, size, length,
                                    probs, pi)
    if size == 1 and simulator:
        freq = _frequencing(tips, site=False)
        ts = tree.format('newick').strip()
        if exp_model:
            s = sim.sim(simulator, ts, model=exp_model, length=length,
                        freq=freq, n=n,  verbose=verbose)

            tree, rates, records, size = _load(s)
            info('Identifying parallel and convergent amino acid replacements in '
                 'simulated sequences.')
            _, par, con, detail = _pc(tree, rates, records, size, length)

            for k, v in par.items():
                pars[k].append(np.mean(v))
                cons[k].append(np.mean(con[k]))
                pars[k].extend(v), cons[k].extend(con[k])
            details.extend(detail)
        else:
            info('No exp_model was assigned, no sequences simulated, expected P&C not calculated.')
        
    if any([pars, cons, details]) and save:
        info('Writing identified parallel and convergent amino acid '
             'replacements to files.')
        bn = utilities.basename(sequence)
        counts = ''.join([bn, '.counts.tsv'])
        changes = ''.join([bn, '.details.tsv'])
        
        with open(counts, 'w') as o, open(changes, 'w') as c:
            o.write('{}\n'.format('\t'.join(h1)))
            s = lambda x: '{:.4f}'.format(x) if isinstance(x, float) else str(x)
            o.writelines('P\t{}\t{}\n'.format(k, '\t'.join([s(x) for x in v]))
                         for k, v in pars.items())
            o.writelines('C\t{}\t{}\n'.format(k, '\t'.join([s(x) for x in v]))
                         for k, v in cons.items())
            
            c.write('{}\n'.format('\t'.join(h2)))
            c.writelines('{}\t{}\t{}\t{}\t{}\t{}\n'.format(*detail)
                         for detail in details)
    
    return pars, cons, details, length


def _tester(obs, exp, values, alpha=0.05):
    """
    One sample T-test to determine whether the observed value is statistically
    significantly different to the expected value.
    
    :param obs: int, observed value.
    :param exps: list or tuple, a sequence of expected values.
    :param alpha: float, significance level.
    :return: float, p value of the T-test.
    """
    
    try:
        from scipy import stats
    except ImportError:
        warn('SciPy package not installed, statistical test aborted.')
        return None
    
    t = stats.ttest_1samp(values, obs)
    return t[1]


def _poisson(obs, exp, values):
    p = poisson.cdf(obs, exp) if obs <= exp else 1 - poisson.cdf(obs, exp)
    return p


def detect(pair=None, pars=None, cons=None, wd='', tester=None, printout=True):
    if pars is None or cons is None:
        wd = wd if wd else os.getcwd()
        if wd and os.path.isdir(wd):
            fn = os.path.join(wd, 'parcon.counts.tsv')
            if os.path.isfile(fn):
                pars, cons = {}, {}
                with open(fn) as f:
                    for line in f:
                        blocks = line.strip().split()
                        if blocks[0] == 'P':
                            pars[blocks[1]] = blocks[2:]
                        elif blocks[0] == 'C':
                            cons[blocks[1]] = blocks[2:]
            else:
                error('Result {} is not a file or does not exist, detect '
                      'aborted.'.format(fn))
                sys.exit(1)
        else:
            error('The wd {} is not a directory or does not exist, detect '
                  'aborted.'.format(wd))
            sys.exit(1)
    elif not isinstance(pars, dict) and not isinstance(cons, dict):
        error('Invalid pars and cons, they need to be dict objects, detect '
              'aborted.')
        sys.exit(1)
    
    if isinstance(pair, str):
        pairs = [pair]
    elif isinstance(pair, (list, tuple)):
        pairs = pair
    else:
        info('No interested branch pair assigned, doing test for all branch pairs.')
        pairs = list(pars.keys())
        
    ps = []
    for p in pairs:
        if p in pars:
            ps.append(p)
        else:
            pr = '-'.join(p.split('-')[::-1])
            if pr in pars:
                ps.append(pr)
            else:
                warn('Invalid branch pair {} was ignored'.format(p))
        
    R = namedtuple('result', 'bp po pe p1 co ce p2')
    results = []
    for item in ps:
        p, c = pars[item], cons[item]
        tester = tester if tester else _poisson
        (po, pe, *pv), (co, ce, *cv) = p, c
        po, pe, pv = int(po), float(pe), [int(i) for i in pv]
        co, ce, cv = int(co), float(ce), [int(i) for i in cv]
        
        p1, p2 = tester(po, pe, pv), tester(co, ce, cv)

        results.append(R(*[item, po, pe, p1, co, ce, p2]))

    if printout and results:
        top = ''.join(['+', '-' * 78, '+'])
        h1 = ''.join(['|', ' ' * 34, '|',
                      'Parallelism'.center(21), '|',
                      'Convergence'.center(21), '|'])
        line1 = ''.join(['+', '-' * 34, '+', '------+', '------+', '-------+',
                        '------+', '------+', '-------+'])
        h2 = ''.join(['|', 'Branch Pair'.center(34), '|']
                     + [' Obs. | Exp. |P-value|'] * 2)
        print('Parallel and convergent amino acid replacements among {} '
              'branch pairs'.format(len(results)))
        print(top)
        print(h1)
        print(line1)
        print(h2)
        for result in results:
            print(line1)
            bp = result.bp.center(34)
            empty = 'None'.center(7)
            po, pe = str(result.po).center(6), '{:5.4f}'.format(result.pe)
            p1 = '{:2.1E}'.format(result.p1) if result.p1 else empty
            co, ce = str(result.co).center(6), '{:5.4f}'.format(result.ce)
            p2 = '{:2.1E}'.format(result.p2) if result.p2 else empty
            print('|{}|'.format('|'.join([bp, po, pe, p1, co, ce, p2])))
        print(top)
    return results


def main():
    des = """identifying parallel and convergent amino acid replacements in
orthologous protein sequences"""
    
    epilog = """
Sequence data file covers a wide range of files and formats:
    
    * sequences: raw protein sequence file, need to be in FASTA format
      and a NEWICK format tree is also required for argument tree.
    * msa: multiple sequence alignment file, need to be in FASTA format
      and a NEWICK format tree is also required for argument tree.
    * ancestors: reconstructed ancestral states file, need to be in tsv
      (tab separated) file, the first line needs to start with #TREE,
      second line need to be a blank line, and the rest lines in the
      file need to be tab separated sequence name (or ID) and amino
      acid sequences.
    * simulations: simulated sequences, need to be in tsv file, the
      first line needs to start with #TREE, second line need to be
      a blank line, each dataset need to be separated by a blank line
      and inside each dataset block, each line should consist of tab
      separated sequence name (or ID) and amino acid sequences.
"""

    formatter = argparse.RawDescriptionHelpFormatter
    parse = argparse.ArgumentParser(description=des, prog='ProtParCon-ProtParCon',
                                    usage='%(prog)s SEQUENCE [OPTIONS]',
                                    formatter_class=formatter, epilog=epilog)
    
    parse.add_argument('SEQUENCE',
                       help='Path to the sequence data file.')
    parse.add_argument('-tree',
                       help='Path to the NEWICK format tree file.')
    parse.add_argument('-aligner',
                       help='Path to the executable of an alignment program')
    parse.add_argument('-ancestor',
                       help='Path to the executable of an ancestral states '
                            'reconstruction program.')
    parse.add_argument('-simulator',
                       help='Path to the executable of an sequence simulation '
                            'program.')
    parse.add_argument('-asr_model', default='JTT',
                       help='Model name or model file for ancestral states '
                            'reconstruction.')
    parse.add_argument('-exp_model',
                       help='Model name or model file for sequence simulation.')
    parse.add_argument('-n', default=100, type=int,
                       help='Number of datasets (or duplicates) should be '
                            'simulated.')

    parse.add_argument('-v', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')
    
    args = parse.parse_args()
    s, tree = args.SEQUENCE, args.tree
    aligner, ancestor = args.aligner, args.ancestor
    simulator = args.simulator
    asr_model, exp_model = args.asr_model, args.exp_model
    n, verbose = args.n, args.v
    
    imc(s, tree=tree, aligner=aligner, ancestor=ancestor, simulator=simulator,
        asr_model=asr_model, exp_model=exp_model, n=n, verbose=verbose,
        save=True)


if __name__ == '__main__':
    main()
