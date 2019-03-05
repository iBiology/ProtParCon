#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Providing a common interface for simulating amino acid sequences using
various simulation programs. At this stage, the module only support simulate
sequences using EVOLVER (inside PAML program) and Seq-Gen.

The minimum requirement of this interface asks users to provide a NEWICK format
tree string or tree file and a executable of a simulation program. Users can
also pass additional options to simulate amino acid sequences under various
cases.

Users are recommended to only use function sim() and avoid to use any
private functions inside the module. However, users are strongly recommended to
implement new private functions for additional simulation programs and wrap
them into the general use function aut().
"""

import os
import re
import sys
import shutil
import random
import logging
import tempfile
import argparse

from io import StringIO
try:
    from textwrap import indent
except ImportError:
    from ProtParCon.utilities import indent
from collections import Counter
from subprocess import PIPE, Popen

from ProtParCon.utilities import basename, modeling, Tree
from ProtParCon.models import models
from Bio import Phylo, AlignIO

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error

MODELS = os.path.join(os.path.dirname(__file__), 'ProtParCon', 'data')
AMINO_ACIDS = 'ARNDCQEGHILKMFPSTWYV'
MC_DAT = """0          * 0: paml format (mc.paml); 1:paup format (mc.nex)
{}         * random number seed (odd number)
{} {} {}   * <# seqs>  <# sites>  <# replicates>
-1         * <tree length, use -1 if tree below has absolute branch lengths>
{}
{} {}      * <alpha>  <#categories for discrete gamma>'
{} {}      * <model> [aa substitution rate file, need only if model=2 or 3]
{}
{}
A R N D C Q E G H I
L K M F P S T W Y V
"""


def _guess(exe):
    """
    Guess the name of a sequence simulation program according to its
    executable.

    :param exe: str, path to the executable of an simulation program.
    :return: tuple, name of the simulation program and the corresponding
        function.
    """
    
    wd = tempfile.mkdtemp()
    try:
        process = Popen([exe, '-h'], cwd=wd, stdout=PIPE, stderr=PIPE,
                        universal_newlines=True)
        try:
            outs, errs = process.communicate(timeout=3)
            out = outs or errs
            if 'seq-gen' in out:
                return 'seq-gen', _seqgen
            else:
                return 'evolver', _evolver
        except Exception as exc:
            process.terminate()
            process.wait(timeout=10)
            return 'evolver', _evolver
    except OSError:
        error("The exe ({}) is empty or may not be an valid executable of a "
              "sequence simulation program.".format(exe))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)


def _evolver_parse(wd):
    """
    Parse the work directory of EVOLVER to get the simulated results.
    
    :param wd: str, work directory of EVOLVER.
    :return: tuple, a list of simulated sequences (dict) and a tree object.
    """
    
    simulations, tree = [], None

    notation, ts, ancs, leaves = '', [], [], []
    with open(os.path.join(wd, "ancestral.txt")) as f:
        for line in f:
            line = line.strip()
            if line and line[0].isdigit() and line[-1].isdigit():
                notation = line
                break
    
        for line in f:
            if line.startswith('['):
                break
    
        for line in f:
            if line.strip() and not line.startswith('['):
                ancs.append(line)

    with open(os.path.join(wd, 'mc.paml')) as f:
        for line in f:
            if line.strip():
                leaves.append(line)

    with open(os.path.join(wd, 'simulation.log')) as f:
        for line in f:
            if line.startswith('(') and line.strip().endswith(';'):
                ts.append(line)

    branches = {}
    for branch in notation.split():
        head, tail = branch.split('..')
        branches[tail] = head

    tree, code = [Phylo.read(StringIO(t.replace(' ', '')), 'newick') for t
                  in ts]
    root = None
    for clade, node in zip(tree.find_clades(order='postorder'),
                           code.find_clades(order='postorder')):
        try:
            parent = tree.get_path(clade)[-2]
            parent.name = branches[node.name] if node.name else branches[
                clade.name]
        except IndexError:
            if len(tree.get_path(clade)) == 1:
                root = branches[node.name] if node.name else branches[
                    clade.name]
            else:
                clade.name = root

    number, maps = tree.count_terminals(), {}
    for clade in tree.find_clades():
        if not clade.is_terminal():
            number += 1
            old, new = clade.name, 'NODE{}'.format(number)
            maps[old] = new
            clade.name = new

    if ancs and leaves:
        for inter, leaf in zip(
                AlignIO.parse(StringIO(''.join(ancs)), 'phylip-relaxed'),
                AlignIO.parse(StringIO(''.join(leaves)), 'phylip-relaxed')):
            leaf.extend(inter)
            for record in leaf:
                name = record.id.replace('node', 'NODE')
                record.id = maps.get(name, name)
            simulations.append(leaf)
    return simulations, tree
        
        
def _evolver(exe, tree, length, freq, model, n, seed, gamma, alpha, invp,
             outfile):
    """
    Sequence simulation via EVOLVER.
    
    :param exe: str, path to the executable of EVOLVER.
    :param tree: str, path to the tree (must has branch lengths and in NEWICK
        format).
    :param length: int, the number of the amino acid sites need to be simulated.
    :param freq: list or None, base frequencies of 20 amino acids.
    :param model: str, name of a model a filename of a model file.
    :param n: int, number of datasets (or duplicates) need to be simulated.
    :param seed: int, the seed used to initiate the random number generator.
    :param gamma: int, 0 means discrete Gamma model not in use, any
        positive integer larger than 1 will invoke discrete Gamma model and
        set the number of categories to gamma.
    :param alpha: float, the Gamma shape parameter alpha, without setting, the
        value will be estimated by the program, in case an initial value is
        needed, the initial value of alpha will be set to 0.5.
    :param invp: float, proportion of invariable site.
    :param outfile: pathname of the output ML tree. If not set, default name
        [basename].[program].ML.newick, where basename is the filename of the
        sequence file without extension, program is the name of the ML inference
        program, and newick is the extension for NEWICK format tree file.
    :return: str, path to the simulation output file.
    """
    
    wd = os.path.dirname(outfile)
    cwd = tempfile.mkdtemp(dir=wd)
    dat = 'MCaa.dat'
    tree = Tree(tree, leave=True)
    tn, ts = tree.leaves, tree.string()
    m = modeling(model)
    if m.type == 'custom':
        mf = m.name
    else:
        name = m.name
        if name.lower() in models:
            info('Using {} model for ancestral states reconstruction.'.format(
                    name))
            with open(os.path.join(cwd, name), 'w') as o:
                o.write(models[name.lower()])
            mf = name
        else:
            error('PAML (evolver) does not support model {}.'.format(name))
            sys.exit(1)

    if freq is None:
        mn = 2
        f1, f2 = '', ''
    else:
        mn = 3
        f1 = ' '.join([str(i) for i in freq[:10]])
        f2 = ' '.join([str(i) for i in freq[10:]])

    with open(os.path.join(cwd, dat), 'w') as o:
        o.write(MC_DAT.format(seed, tn, length, n, ts, alpha, gamma,
                              mn, mf, f1, f2))
        
    try:
        info('Simulating sequences using EVOLVER.')
        log = os.path.join(cwd, 'simulation.log')
        with open(log, 'w') as stdoe:
            process = Popen([exe, '7', dat], cwd=cwd, stdout=stdoe,
                            stderr=stdoe, universal_newlines=True)
            code = process.wait()
        if code:
            with open(log) as handle:
                error('Sequence simulation via EVOLVER failed for {} due to:'
                      '\n{}'.format(tree, indent(handle.read(), prefix='\t')))
                sys.exit(1)
        else:
            info('Parsing and saving simulation results.')
            simulations, tree = _evolver_parse(cwd)
            
            try:
                with open(outfile, 'w') as o:
                    o.write('#TREE\t{}\n'.format(tree.format('newick')))
                    for simulation in simulations:
                        o.writelines('{}\t{}\n'.format(s.id, s.seq)
                                     for s in simulation)
                        o.write('\n')
            except OSError:
                error('Failed to save simulation results to {} ('
                      'IOError, permission denied).'.format(outfile))
                outfile = ''
            info('Successfully saved simulation results to {}'.format(outfile))
    except OSError:
        error('Invalid PAML (EVOLVER) executable {}, running EVOLVER failed '
              'for {}.'.format(exe, tree))
        sys.exit(1)
    finally:
        shutil.rmtree(cwd)
    return outfile


def _seqgen(exe, tree, length, freq, model, n, seed, gamma, alpha, invp,
            outfile):
    """
    Sequence simulation via EVOLVER.

    :param exe: str, path to the executable of EVOLVER.
    :param tree: str, path to the tree (must has branch lengths and in NEWICK
        format).
    :param length: int, the number of the amino acid sites need to be simulated.
    :param freq: list or None, base frequencies of 20 amino acids.
    :param model: str, name of a model a filename of a model file.
    :param n: int, number of datasets (or duplicates) need to be simulated.
    :param seed: int, the seed used to initiate the random number generator.
    :param gamma: int, 0 means discrete Gamma model not in use, any
        positive integer larger than 1 will invoke discrete Gamma model and
        set the number of categories to gamma.
    :param alpha: float, the Gamma shape parameter alpha, without setting, the
        value will be estimated by the program, in case an initial value is
        needed, the initial value of alpha will be set to 0.5.
    :param invp: float, proportion of invariable site.
    :param outfile: pathname of the output ML tree. If not set, default name
        [basename].[program].ML.newick, where basename is the filename of the
        sequence file without extension, program is the name of the ML inference
        program, and newick is the extension for NEWICK format tree file.
    :return: str, path to the simulation output file.
    """
    
    wd = os.path.dirname(outfile)
    cmd = [exe, '-l{}'.format(length), '-n{}'.format(n), '-z{}'.format(seed)]
    m = modeling(model)
    if m.type == 'custom':
        try:
            with open(m.name) as handle:
                lines = handle.readlines()
        except IndexError:
            error('Invalid model file {}, Line 22 (amino acid frequencies)'
                  'does not exist in model file.'.format(m.name))
            sys.exit(1)
        r = [line.strip() for line in lines[:19]]
        r = re.sub('\s+', ',', ','.join(r))
        cmd.append('-r{}'.format(r))
        if freq is None:
            freq = re.sub(r'\s+', ',', ','.join(lines[21]))
            cmd.append('-f{}'.format(freq))
            freq = None
    else:
        cmd.append('-m{}'.format(m.name.upper()))
    if freq:
        cmd.append('-f{}'.format(','.join([str(i) for i in freq])))

    if gamma:
        cmd.append('-g{}'.format(gamma))
    if alpha:
        cmd.append('-a{}'.format(alpha))
    cmd.extend(['-wa', '-q'])
        
    cwd = tempfile.mkdtemp(dir=wd)
    output = os.path.join(cwd, 'output.phylip')
    tree = Tree(tree).file(os.path.join(cwd, 'tree.newick'))

    try:
        info('Simulating sequences using Seq-Gen.')
        stdout, stdin = open(output, 'w'), open(tree)
        process = Popen(cmd, cwd=cwd, stdout=stdout, stdin=stdin,
                        stderr=PIPE, universal_newlines=True)
        code = process.wait()
        stdout.close(), stdin.close()
        if code:
            msg = indent(process.stderr.read(), prefix='\t')
            error('Sequence simulation via Seq-Gen failed due to:'
                  '\n{}'.format(tree, msg))
            sys.exit(1) 
        else:
            info('Parsing and saving simulation results.')
            tree = Phylo.read(tree, 'newick')
            number, nodes = tree.count_terminals(), []
            for clade in tree.find_clades():
                if not clade.is_terminal():
                    number += 1
                    clade.name = 'NODE{}'.format(number)
                    nodes.append(str(number))
            try:
                with open(outfile, 'w') as o:
                    o.write('#TREE\t{}\n'.format(tree.format('newick')))
                    with open(output) as f:
                        for line in f:
                            if line.strip():
                                i, s = line.strip().split()
                                if i.isdigit() and s.isdigit():
                                    o.write('\n')
                                else:
                                    if i.isdigit():
                                        i = 'NODE{}'.format(i)
                                    o.write('{}\t{}\n'.format(i, s))
            except OSError:
                error('Failed to save simulation results to {} ('
                      'IOError, permission denied).'.format(outfile))
                outfile = ''
            info('Successfully saved simulation results to {}'.format(outfile))
    except OSError:
        error('Invalid Seq-Gen executable {}, running Seq-Gen failed for '
              '{}.'.format(exe, tree))
        sys.exit(1)
    finally:
        shutil.rmtree(cwd)
    return outfile

        
def sim(exe, tree, model='JTT', msa='', length=100, freq='empirical', n=100,
        seed=0, gamma=4, alpha=0.5, invp=0, outfile='', verbose=False):
    """
    Sequence simulation via EVOLVER.

    :param exe: str, path to the executable of EVOLVER.
    :param tree: str, path to the tree (must has branch lengths and in NEWICK
        format).
    :param length: int, the number of the amino acid sites need to be simulated.
    :param freq: list or None, base frequencies of 20 amino acids.
    :param model: str, name of a model a filename of a model file.
    :param n: int, number of datasets (or duplicates) need to be simulated.
    :param seed: int, the seed used to initiate the random number generator.
    :param gamma: int, 0 means discrete Gamma model not in use, any
        positive integer larger than 1 will invoke discrete Gamma model and
        set the number of categories to gamma.
    :param alpha: float, the Gamma shape parameter alpha, without setting, the
        value will be estimated by the program, in case an initial value is
        needed, the initial value of alpha will be set to 0.5.
    :param invp: float, proportion of invariable site.
    :param outfile: pathname of the output ML tree. If not set, default name
        [basename].[program].ML.newick, where basename is the filename of the
        sequence file without extension, program is the name of the ML inference
        program, and newick is the extension for NEWICK format tree file.
    :return: str, path to the simulation output file.
    """
    
    logger.setLevel(logging.INFO if verbose else logging.ERROR)
    
    t = Tree(tree, leave=True)
    if not t.length:
        error('Unscaled tree: {}, cannot simulate sequences without branch'
              'lengths.'.format(tree))
        sys.exit(1)
    
    fs = None
    if msa:
        if os.path.isfile(msa):
            aln = AlignIO.read(msa, 'fasta')
            length = aln.get_alignment_length()
            if freq == 'estimate':
                counts = Counter({a: 0 for a in AMINO_ACIDS})
                for a in aln:
                    counts += Counter(a.seq)
                total = float(sum([v for k, v in counts.items()
                                  if k in AMINO_ACIDS]))
                if total:
                    fs = ['{:.8f}'.format(counts.get(aa, 0) / total)
                          for aa in AMINO_ACIDS]
                else:
                    warn('Calculate amino acid frequencies failed, zero amino'
                         'acid found, use empirical frequency instead.')
                    fs = None
        else:
            error('MSA {} is not a file or does not exist.'.format(msa))
            sys.exit(1)
    elif length:
        try:
            length = int(length)
        except ValueError:
            error('Invalid length, length should be a integer.')
            sys.exit(1)
    else:
        error('Neither length or msa was provided, simulation aborted.')
        sys.exit(1)

    if freq:
        if freq == 'equal':
            fs = ['0.05'] * 20
        elif freq == 'estimate':
            if not fs:
                warn('Failed to get observed amino acid frequency, use model '
                     'frequency instead.')
        elif freq == 'empirical':
            fs = None
        elif freq.startswith('0') and freq.count(',') == 19:
            fs = [i.strip() for i in freq.split(',')]
        else:
            warn('Unknown frequency encounter, use model frequency instead.')
            fs = None
            
    try:
        seed = int(seed) if seed else random.randint(0, 10000)
    except ValueError:
        warn('Invalid seed, generating seed using random number generator.')
        seed = random.randint(0, 10000)
    
    name, func = _guess(exe)
    if not outfile:
        outfile = os.path.join(os.getcwd(), '{}.simulations.tsv'.format(name))
    return func(exe, tree, length, fs, model, n, seed, gamma, alpha,
                invp, outfile)
    

def main():
    des = 'Perform topology test (AU TEST) over protein alignments.'
    epilog = """
The minimum requirement for running iMC-aut is a multiple sequence alignment
(MSA) file, a tree file and an executable of a AU TEST program.

If users only provide one tree in the tree file, a maximum-likelihood tree will
be inferred and AU-TEST will be performed to test the difference between them.
If a set of trees in NEWICK format was provided in the tree file, only these
trees will be evaluated without inferring the ML tree. In both case, only
the p-value of AU TEST for the first tree will be printed out.
"""
    formatter = argparse.RawDescriptionHelpFormatter
    parse = argparse.ArgumentParser(description=des, prog='ProtParCon-sim',
                                    usage='%(prog)s EXE TREE [OPTIONS]',
                                    formatter_class=formatter, epilog=epilog)

    parse.add_argument('EXE',
                       help='Pathname of the executable of the sequence '
                            'simulation program.')
    parse.add_argument('TREE',
                       help='Pathname of a tree file or a tree string in '
                            'newick format.')
    parse.add_argument('-msa',
                       help='Pathname of the MSA (sequence) file.')
    parse.add_argument('-len', default=100,
                       help='the number of the amino acid sites need to '
                            'be simulated')
    parse.add_argument('-freq', default='empirical',
                       help='Comma separated state frequencies.')
    parse.add_argument('-model', default='JTT',
                       help='Name of the evolutionary model or filename of the '
                            'model file.')
    parse.add_argument('-gamma', default=4,
                       help='The number of categories for the discrete gamma '
                            'rate heterogeneity model.')
    parse.add_argument('-alpha', default=0.5,
                       help='The shape (alpha) for the gamma rate '
                            'heterogeneity.')
    parse.add_argument('-num', default=3,
                       help='Number of datasets to be simulated.')
    parse.add_argument('-seed',
                       help='Random number seed.')
    parse.add_argument('-invp', default=0,
                       help='The proportion of invariable sites.')
    parse.add_argument('-o',
                       help='Pathname of the output file for storing '
                            'simulated sequences.')
    parse.add_argument('-v', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')
    
    args = parse.parse_args()
    exe, tree, msa, length = args.EXE, args.TREE, args.msa, args.len
    freq, model, gamma, alpha = args.freq, args.model, args.gamma, args.alpha
    num, seed, invp, o, verbose = args.num, args.seed, args.invp, args.o, args.v

    sim(exe, tree, model, msa=msa, length=length, freq=freq, n=num, seed=seed,
        gamma=gamma, alpha=alpha, invp=invp, outfile=o, verbose=verbose)


if __name__ == '__main__':
    main()
