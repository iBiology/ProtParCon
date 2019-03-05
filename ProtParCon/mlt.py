#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Providing a common interface for inferring phylogenetic trees from alignments
of protein sequences using maximum-likelihood method.

The minimum requirement for using this module is a multiple sequence alignment
(MSA) file and an executable of a maximum-likelihood tree inference program.
If user wants to take full control of the tree inference, optional parameters
are also accepted depend on the program. Without providing an evolutionary
model, the default model (JTT model for FastTree and LG model for PhyML) or a
model decided by model selection procedure (RAxML and IQ-TREE) will be used.

Users are recommended to only use function ``mlt()`` and avoid to use any
private functions in this module. However, users are strongly recommended to
implement new private functions for additional maximum-likelihood tree
inference program that they are interested and incorporate them into the
general use function ``mlt()``.

.. TODO::
    It seems FastTree does not accept FASTA format file which has a space
    between '>' and the sequence name (or ID), should check FASTA file before
    feed in.
"""

import os
import re
import sys
import random
import shutil
import logging
import tempfile
import argparse
try:
    from textwrap import indent
except ImportError:
    from utilities import indent
from collections import namedtuple
from subprocess import PIPE, Popen, DEVNULL

from Bio import Phylo, AlignIO
from utilities import basename, modeling

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error

FASTTREE_MODELS = ['jtt', 'lg', 'wag']
IQTREE_MODELS = ['BLOSUM62', 'cpREV', 'Dayhoff', 'DCMut', 'FLU', 'HIVb', 'HIVw',
                 'JTT', 'JTTDCMut', 'LG', 'mtART', 'mtMAM', 'mtREV', 'mtZOA',
                 'mtMet', 'mtVer', 'mtInv', 'Poisson', 'PMB', 'rtREV', 'VT',
                 'WAG', 'GTR20']
RAXML_MODELS = ['DAYHOFF', 'DCMUT', 'JTT', 'MTREV', 'WAG', 'RTREV', 'CPREV',
                'VT', 'BLOSUM62', 'MTMAM', 'LG', 'MTART', 'MTZOA', 'PMB',
                'HIVB', 'HIVW', 'JTTD', 'CMUT', 'FLU', 'STMTREV', 'DUMMY',
                'DUMMY2', 'AUTO', 'LG4M', 'LG4X', 'PROT_FILE', 'GTR_UNLINKED',
                'GTR']
PHYML_MODELS = ['Blosum62', 'CpREV', 'DCMut', 'Dayhoff', 'HIVb', 'HIVw', 'JTT',
                'MtArt', 'MtMam', 'MtREV', 'RtREV', 'VT', 'WAG']

MODEL = namedtuple('model', 'name frequency gamma rates invp type')


def _guess(exe):
    """
    Guess the name of a tree inference program according to its executable.
    
    :param exe: str, path to the executable of an maximum-likelihood tree
        inference program.
    :return: tuple, name of the program and the corresponding function.
    """
    
    programs, func = {'FastTree': _fasttree, 'IQ-TREE': _iqtree,
                      'RAxML': _raxml, 'PhyML': _phyml}, None
    if exe:
        try:
            process = Popen([exe, '-h'], stdout=PIPE, stderr=PIPE,
                            universal_newlines=True)
            outs, errs = process.communicate(timeout=10)
            out = outs or errs
            for name in programs.keys():
                if name in out[:50]:
                    return name, programs[name]
            error('Unknown exe: {}, the exe may not be an acceptable tree '
                  'inference program, exited.')
            sys.exit(1)
        except OSError:
            error('Invalid exe: {}, exe does not exist, exited.'.format(exe))
            sys.exit(1)
    else:
        error('The exe of a tree inference program is empty, exited.')
        sys.exit(1)


def _fasttree(exe, msa, model, cat, gamma, alpha, freq, invp, start_tree,
              constraint_tree, seed, outfile):
    """
    Infer ML phylogenetic tree using FastTree.
    """
    
    if model.type == 'custom':
        error('Invalid model, FastTree does not accept custom model file '
              '{}.'.format(model))
        sys.exit(1)
    else:
        name = model.name if model.name else 'JTT'
        if name.lower() == 'jtt':
            info('Inferring ML tree using JTT model.')
            m = ''
        else:
            if name.lower in ['lg', 'wag']:
                info('Inferring ML tree using {} model.'.format(
                        model.upper()))
                m = '-{}'.format(name)
            else:
                error('Invalid model {}, FastTree only accepts JTT, WAG, '
                      'and LG model.'.format(name))
                sys.exit(1)
    
    info('Inferring ML tree for {} using FastTree.'.format(msa))
    args = [exe, '-quiet', '-nopr', '-seed', str(seed)]
    if start_tree and os.path.isfile(start_tree):
        args.extend(['-intree', start_tree])
    if m:
        args.append(m)
    if gamma or model.gamma:
        args.append('-gamma')
    if cat:
        if cat != 20:
            args.extend(['-cat', str(cat)])
    else:
        args.append('-nocat')
    args.append(os.path.basename(msa))
    
    tree = outfile if outfile else '{}.FastTree.ML.newick'.format(basename(msa))
    try:
        info('Running FastTree using the following command:\n\t'
             '{}'.format(' '.join(args)))
        with open(tree, 'w') as stdout:
            process = Popen(args, cwd=os.path.dirname(msa), stdout=stdout, stderr=PIPE,
                            universal_newlines=True)
        code = process.wait()
        if code:
            msg = indent(process.stderr.read(), prefix='\t')
            print(msg)
            error('Inferring ML tree failed for {}\n{}'.format(msa, msg))
            if os.path.isfile(tree):
                os.remove(tree)
            sys.exit(1)
        else:
            info('Successfully save inferred ML tree to {}.'.format(tree))
    except OSError:
        error('Inferring ML tree failed for {}, executable (exe) {} of '
              'FastTree is empty or invalid.'.format(msa, exe))
        if os.path.isfile(tree):
            os.remove(tree)
        sys.exit(1)
    return tree


def _iqtree(exe, msa, model, cat, gamma, alpha, freq, invp, start_tree,
            constraint_tree, seed, outfile):
    """
    Infer ML phylogenetic tree using IQ-TREE.
    """
    
    wd = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(msa)))
    shutil.copy(msa, os.path.join(wd, 'msa.fa'))
    if model.name:
        name = model.name.upper()
        if model.type == 'builtin':
            frequency = freq or model.frequency
            if frequency == 'estimate':
                frequency = 'FO'
            else:
                frequency = 'F'
            
            invpro = 'I' if (invp or model.invp) else ''
            gamma = gamma or model.gamma
            gamma = 'G{}'.format(gamma) if gamma else ''
            cat = cat if cat not in ('None', 'none', None) else 0
            rates = cat or model.rates
            rates = 'R{}'.format(rates) if rates else ''
            m = ['-m', '+'.join((i for i in [name, frequency, invpro, gamma, rates] if i))]
        else:
            m = ['-mdef', model.name]
    else:
        m = ['-m', 'TEST']

    info('Inferring ML tree for {} using IQ-TREE.'.format(msa))
    # iqtree_cmd = 'exe -s seq -t user_tree -pre prefix -nt 1 -seed seed
    # -quiet -m MFP -g constraint_tree'
    args = [exe, '-s', 'msa.fa', '-nt', '1', '-pre', 'tmf', '-seed', str(seed)]
    args.extend(m)
    # args.append('-quiet')
    
    if gamma and alpha:
        args.extend(['-a', str(alpha)])
    if invp:
        args.extend(['-i', str(invp)])
    if constraint_tree:
        args.extend(['-g', constraint_tree])
    try:
        # info('Running IQTree using the following command:\n\t'
        #      '{}'.format(' '.join(args)))
        process = Popen(args, cwd=wd, stdout=PIPE, stderr=PIPE,
                        universal_newlines=True)
        code = process.wait()
        if code:
            msg = indent(process.stderr.read(), prefix='\t')
            error('Inferring ML tree failed for {}\n{}'.format(msa, msg))
            sys.exit(1)
        else:
            tree = outfile if outfile else '{}.IQ-TREE.ML.newick'.format(
                    basename(msa))
            try:
                tree = shutil.copy(os.path.join(wd, 'tmf.treefile'), tree)
                info('Successfully save inferred gene tree to {}.'.format(
                        tree))
            except OSError:
                error('Path of outfile {} is not writeable, saving tree to '
                      'file failed.'.format(tree))
                sys.exit(1)
    
    except OSError:
        error('Inferring tree failed for {}, executable (exe) {} of '
              'IQ-TREE is invalid.'.format(msa, exe))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)
    
    return tree


def _raxml(exe, msa, model, cat, gamma, alpha, freq, invp, start_tree,
           constraint_tree, seed, outfile):
    """
    Infer ML phylogenetic tree using RAxML.
    """

    wd = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(msa)))
    shutil.copy(msa, os.path.join(wd, 'msa.fa'))
    if model.type == 'builtin':
        frequency = freq or model.frequency
        frequency = 'X' if frequency == 'estimate' else 'F'
        gamma = gamma or model.gamma
        invpro = 'I' if (model.invp or invp) else ''
        gamma = 'GAMMA' if gamma else 'CAT'
        if model.name:
            mm = ''.join(['PROT', gamma, invpro, model.name.upper(), frequency])
            model = ['-m', mm]
        else:
            model = ['-m', ''.join(['PROT', invpro, gamma, 'AUTO'])]
    else:
        model = ['-P', model.name]
    
    info('Inferring ML tree for {} using RAxML.'.format(msa))
    args = [exe, '-s', 'msa.fa', '-n', 'RAxML', '-p', str(seed), '--silent']
    args.extend(model)
    cat = cat if cat not in ('None', 'none', None) else 0
    if cat:
        args.extend(['-c', str(cat)])
    if start_tree:
        args.extend(['-t', start_tree])
    if constraint_tree:
        args.extend(['-r', constraint_tree])
    try:
        # info('Running FastTree using the following command:\n\t'
        #      '{}'.format(' '.join(args)))
        process = Popen(args, cwd=wd, stdout=DEVNULL, stderr=PIPE,
                        universal_newlines=True)
        code = process.wait()
        if code:
            msg = indent(process.stderr.read(), prefix='\t')
            error('Tree inferring failed for {}\n{}'.format(msa, msg))
            tree = ''
        else:
            tree = outfile if outfile else '{}.RAxML.ML.newick'.format(
                    basename(msa))
            try:
                tree = shutil.copy(os.path.join(wd, 'RAxML_bestTree.RAxML'),
                                   tree)
                info('Successfully save inferred ML tree to {}.'.format(
                        tree))
            except OSError:
                error('Path of outfile {} is not writeable, saving tree to '
                      'file failed.'.format(tree))
                sys.exit(1)
    except OSError:
        error('Tree inferring failed for {}, executable (exe) {} is '
              'invalid.'.format(msa, exe))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)
    
    return tree


def _phyml(exe, msa, model, cat, gamma, alpha, freq, invp, start_tree,
           constraint_tree, seed, outfile):
    """
    Infer ML phylogenetic tree using PhyML.
    """
    
    # cmd = 'exe -i seq -d aa -m JTT -f e|m -v invariable -c 4 -a gamma-alpha
    # --quiet --r_seed num -u user_tree_file'
    wd = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(msa)))
    alignment = 'temporary.alignment.phylip'
    AlignIO.convert(msa, 'fasta', os.path.join(wd, alignment), 'phylip')
    
    if model.type == 'builtin':
        m = ['-m', model.name] if model.name else ['-m', 'LG']
    else:
        m = ['-m', 'custom', '--aa_rate_file', model.name]
    
    info('Inferring ML tree for {} using PhyML.'.format(msa))
    args = [exe, '-i', alignment, '-d', 'aa', '--r_seed', str(seed), '--quiet']
    args.extend(m)
    cat = cat if cat not in ('None', 'none', None) else 0
    if cat:
        args.extend(['--free_rates', cat])
    gamma = gamma or model.gamma
    if gamma:
        args.extend(['-c', str(gamma)])
    if alpha:
        args.extend(['-a', str(alpha)])
    frequency = freq or model.frequency
    frequency = 'X' if frequency == 'estimate' else 'F'
    if frequency == 'estimate':
        args.extend(['-f', 'e'])
    else:
        args.extend(['-f', 'm'])
    if start_tree:
        args.extend(['-u', start_tree])
        if constraint_tree:
            args.extend(['--constraint_file', constraint_tree])
    if model.invp:
        if invp:
            args.extend(['-v', str(invp)])
        else:
            args.extend(['-v', 'e'])
    else:
        if invp:
            args.extend(['-v', str(invp)])
    
    try:
        # info('Running FastTree using the following command:\n\t'
        #      '{}'.format(' '.join(args)))
        process = Popen(args, cwd=wd, stdout=PIPE, stderr=PIPE,
                        universal_newlines=True)
        code = process.wait()
        if code:
            msg = process.stderr.read() or process.stdout.read()
            msg = indent(msg, prefix='\t')
            error('Tree inferring failed for {}\n{}'.format(msa, msg))
            sys.exit(1)
        else:
            tree = outfile if outfile else '{}.PhyML.ML.newick'.format(
                    basename(msa))
            try:
                out = '{}{}'.format(os.path.join(wd, alignment),
                                    '_phyml_tree.txt')
                tree = shutil.copy(out, tree)
                info('Successfully save inferred ML tree to {}.'.format(
                        tree))
            except OSError:
                error('Path of outfile {} is not writeable, saving tree to '
                      'file failed.'.format(tree))
                sys.exit(1)
    except OSError:
        error('Tree inferring failed for {}, executable (exe) {} of PhyML '
              'is invalid.'.format(msa, exe))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)
    return tree


def mlt(exe, msa, model='', cat=0, gamma=0, alpha=0.0, freq='empirical',
        invp=0.0, start_tree='', constraint_tree='', seed=0, outfile='',
        verbose=False):
    
    """
    Common interface for inferring ML phylogenetic tree.

    :param msa: str, path of the multiple sequence alignment (FASTA) file.
    :param exe: str, path of the executable of the ML tree inference program.
    :param model: str, name of the model or path of the model file.
    :param cat: int, invoke rate heterogeneity (CAT) model and set the
        number of categories to the corresponding number cat, if CAT model is
        not in use, it will be ignored. When FastTree is in use, set cat=None
        to invoke nocat mode.
    :param gamma: int, 0 means discrete Gamma model not in use, any
        positive integer larger than 1 will invoke discrete Gamma model and
        set the number of categories to gamma.
    :param alpha: float, the Gamma shape parameter alpha, without setting, the
        value will be estimated by the program, in case an initial value is
        needed, the initial value of alpha will be set to 0.5.
    :param freq: str, the base frequencies of the twenty amino acids.
        Accept empirical, or estimate, where empirical will set frequencies
        use the empirical values associated with the specified substitution
        model, and estimate will use a ML estimate of base frequencies.
    :param invp: float, proportion of invariable site.
    :param start_tree: str, path of the starting tree file, the tree file
        must be in NEWICK format.
    :param constraint_tree: str, path of the constraint tree file, the tree
        file muse be in NEWICK format.
    :param seed: int, the seed used to initiate the random number generator.
    :param outfile: pathname of the output ML tree. If not set, default name
        [basename].[program].ML.newick, where basename is the filename of the
        sequence file without extension, program is the name of the ML inference
        program, and newick is the extension for NEWICK format tree file.
    :param verbose: bool, invoke verbose or silent process mode, default:
        False, silent mode.
    :return: path of the maximum-likelihood tree file.
    """
    
    level = logging.INFO if verbose else logging.ERROR
    logger.setLevel(level)

    if not os.path.isfile(msa):
        error('Alignment {} is not a file or does not exist.'.format(msa))
        sys.exit(1)
    
    program, func = _guess(exe)

    if program == 'FastTree':
        if cat in ('None', 'none', None):
            cat = 0
        else:
            cat = 20
    
    try:
        seed = int(seed) if seed else random.randint(0, 10000)
    except ValueError:
        warn('Invalid seed, generating seed using random number generator.')
        seed = random.randint(0, 10000)
    
    model = modeling(model)
    tree = func(exe, msa, model=model, cat=cat, gamma=gamma, alpha=alpha,
                freq=freq, invp=invp, start_tree=start_tree, seed=seed,
                constraint_tree=constraint_tree, outfile=outfile)
    return tree


def main():
    des = 'Infer maximum likelihood phylogenetic tree from protein alignments.'
    epilog = """
The minimum requirement for running iMC-mlt is a multiple sequence alignment
(MSA) file and an executable of a maximum-likelihood tree inference program.
If user want to take full control of the tree inference, optional parameters
are also accepted depend on the program.

Without provide an evolutionary model, the default model (JTT model for FastTree
and LG model for PhyML) or a model decided by model selection procedure (RAxML
and IQ-TREE) will be used.

The model needs to be in the format of MODEL+<FreqType>+<RateType> or a model
(text) file, where MODEL is a model name (e.g. JTT, WAG, ...), FreqType decides
how the model will handle amino acid frequency (e.g. F, empirical AA frequencies
from the data; FO, ML optimized AA frequencies from the data or FQ, Equal AA
frequencies).

Some options (e.g. start tree, constraint tree) may not supported by all of the
four programs, they will be automatically ignored.
"""
    formatter = argparse.RawDescriptionHelpFormatter
    parse = argparse.ArgumentParser(description=des, prog='ProtParCon-mlt',
                                    usage='%(prog)s EXECUTABLE MSA [OPTIONS]',
                                    formatter_class=formatter,
                                    epilog=epilog)
    
    parse.add_argument('EXECUTABLE',
                       help='Pathname of the executable of the tree inference '
                            'program.')
    parse.add_argument('MSA',
                       help='Pathname of the alignment file in fasta format.')
    parse.add_argument('-m', '--model',
                       help='Name of the evolutionary model or filename of the '
                            'model file.')
    parse.add_argument('-c', '--category',
                       help='Invoke CAT model and set the number of '
                            'categories to value of cat.', default=0)
    parse.add_argument('-g', '--gamma', default=0,
                       help='Invoke discrete Gamma model and set the number of '
                            'categories to value of gamma.')
    parse.add_argument('-a', '--alpha', help='The Gamma shape parameter alpha.')
    parse.add_argument('-f', '--frequency',
                       help='Amino acid frequency, either m or e.')
    parse.add_argument('-i', '--invp', help='Proportion of invariable sites.')
    parse.add_argument('-p', '--stree', help='Path of the starting tree file.')
    parse.add_argument('-q', '--ctree',
                       help='Path of the constraint tree file.')
    parse.add_argument('-s', '--seed',
                       help='The seed for initiating the random number '
                            'generator.')
    parse.add_argument('-o', '--output',
                       help='Pathname of the ML tree output file.')
    parse.add_argument('-v', '--verbose', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')

    args = parse.parse_args()
    msa, exe = args.MSA, args.EXECUTABLE

    mlt(exe, msa, model=args.model, cat=args.category, gamma=args.gamma,
        alpha=args.alpha, freq=args.frequency, invp=args.invp,
        start_tree=args.stree, constraint_tree=args.ctree,
        seed=args.seed, verbose=args.verbose)
    

if __name__ == '__main__':
    main()
