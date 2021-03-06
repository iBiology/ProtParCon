#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Providing a common interface for identifying parallel and convergent amino
acid replacements in orthologous protein sequences. In order to make this
module for general use, function ``ProtParCon()`` is built on top of other
modules to facilitate the identification of parallel and convergent amino
acid replacements using a wide range of sequence data. Depending on the
sequence data, optional parameters and external programs may be required.
"""

import os
import sys
import glob
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

from ProtParCon import msa, asr, aut, sim, detect, utilities
from ProtParCon.models import models

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error

AMINO_ACIDS = 'ARNDCQEGHILKMFPSTWYV'


def _pairing(tree, indpairs=True):
    """
    Checking whether two branches are sister branch pair or branch pair sharing
    the same evolutionary path.
    
    :param tree: object, a tree object.
    :param indpairs: bool, only return independent branch pairs if true,
        or return all branch pairs if False.
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
    if indpairs:
        pairs = [(b1, b2) for (b1, b2) in combinations(branches, 2)
                 if comparable(tree, b1, b2)]
    else:
        pairs = [(b1, b2) for (b1, b2) in combinations(branches, 2)]
    return branches, pairs
    

def _load(tsv):
    """
    Load tree, rates, and data blocks from a tsv file.
    
    :param tsv: str, path to the tsv file stores ancestral states or simulated
        sequences.
    :return: tuple, tree, rates (list) and sequence records (defaultdict).
    """

    tree, rates, records, aps = None, [], defaultdict(list), {}
    with open(tsv) as handle:
        for line in handle:
            blocks = line.strip().split()
            if len(blocks) >= 2:
                if blocks[0] == '#TREE' and blocks[1].endswith(';'):
                    tree = Phylo.read(StringIO(blocks[1]), 'newick')
                elif blocks[0] == '#RATES':
                    rates = [float(i) for i in blocks[1:]]
                elif blocks[0].startswith('#NODE'):
                    k = blocks[0].replace('#', '')
                    ps = blocks[1].split(')')[:-1]
                    aps[k] = [p.split('(') for p in ps]
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
    return tree, rates, records, aps, size


def _sequencing(sequence, tree, aligner, ancestor, wd, asr_model, verbose):
    """
    Identify the type of the sequence file.
    
    :param sequence: str, path to a sequence data file.
    :param tree: str, path to a NEWICK tree file.
    :return: tuple, sequence, alignment, ancestor, and simulation data file.
    """
    
    if tree:
        utilities.Tree(tree, leave=True)
        AA, lengths, aa = set(AMINO_ACIDS), [], []

        with open(sequence) as handle:
            line = handle.readline().strip()
            if line.startswith('>'):
                handle.seek(0)
                records = SeqIO.parse(handle, 'fasta')
                for record in records:
                    lengths.append(len(record.seq))
                    aa.append(set(record.seq).issubset(AA))
            else:
                error('NEWICK format tree was provided, but the sequence file '
                      'was not in the FASTA format.')
                sys.exit(1)

        if len(set(lengths)) == 1:
            alignment = sequence
            if all(aa):
                trimmed = alignment
            else:
                trimmed = ''.join([utilities.basename(alignment),
                                   '.trimmed.fasta'])
                if os.path.isfile(trimmed):
                    info('Using pre-existed trimmed alignment file.')
                else:
                    _, trimmed = utilities.trim(alignment, outfile=trimmed)
        else:
            if aligner:
                aler, _ = msa._guess(aligner)
                outfile = ''.join([utilities.basename(sequence),
                                  '.{}.fasta'.format(aler)])
                if os.path.isfile(outfile):
                    info('Using pre-existed alignment file')
                    alignment = outfile
                    trimmed = ''.join(
                            [utilities.basename(alignment), '.trimmed.fasta'])
                    if os.path.isfile(trimmed):
                        info('Using pre-existed trimmed alignment file.')
                    else:
                        _, trimmed = utilities.trim(alignment, outfile=trimmed)
                else:
                    trimmed = msa.msa(aligner, sequence, verbose=verbose,
                                        outfile=outfile, trimming=True)
            else:
                error('FASTA format sequence file was provided, but no '
                      'alignment program was provided.')
                sys.exit(1)
        
        if trimmed:
            if ancestor:
                if trimmed.endswith('.trimmed.fasta'):
                    name = trimmed.replace('.trimmed.fasta', '')
                else:
                    name = trimmed
                
                aser, _ = asr._guess(ancestor)
                outfile = '{}.{}.tsv'.format(utilities.basename(name), aser)
                if os.path.isfile(outfile):
                    info('Using pre-existed ancestral states sequence file.')
                    sequence = outfile
                else:
                    sequence = asr.asr(ancestor, trimmed, tree, asr_model,
                                       verbose=verbose, outfile=outfile)
            else:
                error('No ancestral reconstruction program was provided.')
                sys.exit(1)
        else:
            sys.exit(1)
    
    tree, rate, records, aps, size = _load(sequence)
    return tree, rate, records, aps, size, sequence


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


def _pc(tree, rates, records, aps, size, length, probs, pi, indpairs,
        threshold):
    branches, pairs = _pairing(tree, indpairs=indpairs)
    pars, cons, divs = defaultdict(list), defaultdict(list), defaultdict(list)
    details = []

    detail = namedtuple('replacement', 'category position pair r1 r2 dataset')
    for i in range(size):
        record = {k: v[i] for k, v in records.items()}
        for pair in pairs:
            (t1p, t1), (t2p, t2) = pair
            name = '-'.join([t1, t2])
            po, pe, co, ce, do = 0, 0, 0, 0, 0
            for pos in range(length):
                if threshold and aps:
                    s1p, p_s1p = aps[t1p][pos] if t1p in aps else (record[t1p][
                        pos], 1.0)
                    s1, p_s1 = aps[t1][pos] if t1 in aps else (record[t1][
                        pos], 1.0)
                    s2p, p_s2p = aps[t2p][pos] if t2p in aps else (record[t2p][
                        pos], 1.0)
                    s2, p_s2 = aps[t2][pos] if t2 in aps else (record[t2][
                        pos], 1.0)
                    if not all([True if float(p) >= threshold else False
                                for p in [p_s1p, p_s1, p_s2p, p_s2]]):
                        continue
                else:
                    s1p, s1 = record[t1p][pos], record[t1][pos]
                    s2p, s2 = record[t2p][pos], record[t2][pos]

                if s1p != s1 and s2p != s2:
                    if s1 == s2:
                        if size == 1 and i == 0:
                            label = 'OBSERVATION'
                        else:
                            label = 'SIMULATION-{:05d}'.format(i + 1)
                        r1, r2 = '{}{}'.format(s1p, s1), '{}{}'.format(s2p, s2)
                        if s1p == s2p:
                            po += 1
                            cat = 'P'
                        else:
                            co += 1
                            cat = 'C'
                        details.append(detail(cat, pos, name, r1, r2, label))
                    else:
                        if size == 1 and i == 0:
                            label = 'OBSERVATION'
                        else:
                            label = 'SIMULATION-{:05d}'.format(i + 1)
                        r1, r2 = '{}{}'.format(s1p, s1), '{}{}'.format(s2p, s2)
                        do += 1
                        cat = 'D'
                        details.append(detail(cat, pos, name, r1, r2, label))
                            
                if i == 0 and size == 1 and rates and probs is not None:
                    p, c = _prob(tree, rates, record, pos, pair, probs, pi)
                    pe += p
                    ce += c
            if rates and probs is not None:
                pars[name].extend([po, pe])
                cons[name].extend([co, ce])
                divs[name].extend([do, 0.0])
            else:
                pars[name].append(po)
                cons[name].append(co)
                divs[name].append(do)
    
    return tree, pars, cons, divs, details
    

def _load_matrix(model):
    probs, pi = np.zeros((20, 20)), np.zeros((20,))
    model = utilities.modeling(model).name
    if model.lower() in ('jtt', 'jones'):
        handle = StringIO(models['jtt'])
        model = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                             'ProtParCon', 'data', 'jtt')
    else:
        if os.path.isfile(model):
            handle = open(model)
        else:
            error('Unsupported model for computing expected changes, '
                  'calculation aborted.')
            return None, None
    
    for line in handle:
        fields = line.strip().split()
        if len(fields) == 20 and all([i.replace('.', '').isdigit()
                                      for i in fields]):
            pi = np.array([float(i) for i in fields])
            break

    n = 0
    for line in handle:
        fields = line.strip().split()
        if len(fields) == 20 and all([i.replace('.', '').isdigit()
                                      for i in fields]):
            probs[n, :] = [float(field) for field in fields]
            n += 1
    handle.close()

    return probs, pi
    
    
def imc(sequence, tree='', aligner='', ancestor='', simulator='',
        asr_model='JTT', exp_model='JTT', n=100, divergent=True, indpairs=True,
        threshold=0.0, exp_prob=False, verbose=False):
    
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
    :param asr_model: str, model name or model file for ancestral states
        reconstruction, default: JTT.
    :param exp_model: str, model name or model file for estimate expected
        changes based on simulation or replacement probability manipulation,
        default: JTT.
    :param n: int, number of datasets (or duplicates) should be simulated.
    :param divergent: bool, identify divergent changes if True, or only
        identify parallel and convergent changes if False.
    :param indpairs: bool, only identify changes for independent branch pairs
        if true, or identify changes for all branch pairs if False.
    :param threshold: float, a probability threshold that ranges from 0.0 to
        1.0. If provided, only ancestral states with probability equal or
        larger than the threshold will be used, default: 0.0.
    :param exp_prob: bool, calculate the probability of expected changes if set
        to True and the exp_model contains a probability matrix. Time consuming
        process, be patient for the calculation.
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

    basename = utilities.basename(sequence)
    rs = _sequencing(sequence, tree, aligner, ancestor, wd, asr_model, verbose)
    tree, rates, records, aps, size, sequence = rs

    basename_more = utilities.basename(sequence)
    pars, cons, divs, details, aup = None, None, None, None, None
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
                if exp_prob:
                    probs, pi = _load_matrix(exp_model)
                    if probs is not None:
                        h1.append('EXP')
    else:
        h1.append('EXP')
        h1.extend(['SIM-{}'.format(i + 1) for i in range(size)])

    tips = [v[0] for k, v in records.items() if not k.startswith('NODE')]
    length = len(tips[0])
    if size > 1:
        info('Estimating expected changes ... ')
    else:
        info('Identifying observed changes ...')
    tree, pars, cons, divs, details = _pc(tree, rates, records, aps, size,
                                          length, probs, pi, indpairs,
                                          threshold)

    if size == 1 and simulator:
        freq = _frequencing(tips, site=False)
        ts = tree.format('newick').strip()
        out = '{}.{}.tsv'.format(basename, sim._guess(simulator)[0])

        s = sim.sim(simulator, ts, model=exp_model, length=length,
                    freq=freq, n=n,  outfile=out, verbose=verbose)
        
        if s and os.path.isfile(s):
            tree, rates, records, aps, size = _load(s)
            info('Estimating expected changes ... ')
            tree, par, con, div, detail = _pc(tree, rates, records, aps,
                                                  size, length, None, None,
                                                  indpairs, threshold)

            for k, v in par.items():
                pars[k].append(np.mean(v))
                cons[k].append(np.mean(con[k]))
                divs[k].append(np.mean(div[k]))
                pars[k].extend(v), cons[k].extend(con[k])
                divs[k].extend(div[k])
            details.extend(detail)

    if any([pars, cons, divs, details]):
        info('Writing identified parallel and convergent amino acid '
             'replacements to files.')
        counts = ''.join([basename_more, '.counts.tsv'])
        changes = ''.join([basename_more, '.details.tsv'])
        
        with open(counts, 'w') as o, open(changes, 'w') as c:
            o.write('{}\n'.format('\t'.join(h1)))
            s = lambda x: '{:.4f}'.format(x) if isinstance(x, float) else str(x)
            o.writelines('P\t{}\t{}\n'.format(k, '\t'.join([s(x) for x in v]))
                         for k, v in pars.items())
            o.writelines('C\t{}\t{}\n'.format(k, '\t'.join([s(x) for x in v]))
                         for k, v in cons.items())
            o.writelines('D\t{}\t{}\n'.format(k, '\t'.join([s(x) for x in v]))
                         for k, v in divs.items())
            
            c.write('{}\n'.format('\t'.join(h2)))
            c.writelines('{}\t{}\t{}\t{}\t{}\t{}\n'.format(*detail)
                         for detail in details)
    
    return pars, cons, divs, details, length


def main():
    des = """Identifying parallel and convergent amino acid replacements in
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
    parse = argparse.ArgumentParser(description=des,
                                    prog='imc',
                                    usage='%(prog)s SEQUENCE [OPTIONS]',
                                    formatter_class=formatter, epilog=epilog)
    
    parse.add_argument('SEQUENCE',
                       help='Path to the sequence data file.')
    parse.add_argument('-t', '--tree',
                       help='Path to the NEWICK format tree file.')
    parse.add_argument('-l', '--aligner',
                       help='Path to the executable of an alignment program')
    parse.add_argument('-a', '--ancestor',
                       help='Path to the executable of an ancestral states '
                            'reconstruction program.')
    parse.add_argument('-s', '--simulator',
                       help='Path to the executable of an sequence simulation '
                            'program.')
    parse.add_argument('-m', '--asr_model', default='JTT',
                       help='Model name or model file for ancestral states '
                            'reconstruction.')
    parse.add_argument('-r', '--exp_model', default='JTT',
                       help='Model name or model file for sequence simulation.')
    parse.add_argument('-n', '--number', default=100, type=int,
                       help='Number of datasets (or duplicates) should be '
                            'simulated.')
    parse.add_argument('-p', '--probability', default=0.0, type=float,
                       help='a probability threshold that ranges from 0.0 to '
                            '1.0. If provided, only ancestral states with '
                            'probability equal or larger than the threshold '
                            'will be used, default: 0.0')
    parse.add_argument('-i', '--indpairs', action='store_false',
                       help='Identify changes for all branch pairs.')
    parse.add_argument('-c', '--exp_prob', action='store_true',
                       help='Calculate the probability of expected changes if '
                            'the exp_model contains a probability matrix. '
                            'Highly time consuming, be patient.')
    parse.add_argument('-v', '--verbose', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')
    
    args = parse.parse_args()
    s, tree = args.SEQUENCE, args.tree
    
    imc(s, tree=tree, aligner=args.aligner, ancestor=args.ancestor,
        simulator=args.simulator, asr_model=args.asr_model,
        exp_model=args.exp_model, n=args.number, exp_prob=args.exp_prob,
        threshold=args.probability, indpairs=args.indpairs,
        verbose=args.verbose)


if __name__ == '__main__':
    main()
