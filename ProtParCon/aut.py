#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Providing a common interface for performing topology test (`AU test`_) over
protein alignments using various programs.

Users are asked to provide a multiple sequence alignment (MSA) file, a NEWICK 
format tree file, and a topology test program's executable. If only one tree 
in the tree file was provided, a maximum-likelihood (ML) tree would be inferred 
and AU test will be performed to test the difference between the user specified 
tree and the ML tree. If a set of trees in NEWICK format was provided in the 
tree file, only these trees would be evaluated without reconstructing the ML 
tree. In both cases, only the p-value of AU test for the first tree will be 
returned.

Users are recommended only to use function ``aut()`` and avoid to use any
private functions inside the module. However, users are recommended to
implement new private functions for additional topology test programs that
they are interested and incorporate them into the general use function 
``aut()``.

.. _`AU test`: https://academic.oup.com/sysbio/article/51/3/492/1616895

"""

import os
import sys
import shutil
import random
import logging
import tempfile
import argparse

from io import StringIO
from subprocess import PIPE, Popen
try:
    from textwrap import indent
except ImportError:
    from ProtParCon.utilities import indent

from ProtParCon.utilities import basename, Tree
from ProtParCon.mlt import mlt
from Bio import Phylo

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error


def _iqtree(exe, msa, tree, model, seed):
    """
    Perform topology test (AU test) using IQ-TREE_.
    
    :param msa: str, path to a FASTA format MSA file.
    :param tree: str, path to a NEWICK format tree file.
    :param exe: str, path to a topology test program's executable.
    :param model: str, name of the model or path of the model file.
    :param seed: int, the seed used to initiate the random number generator.
    :return: tuple, p-value of the AU test (float), first tree (string), and
        second tree (string).
    
    .. _IQ-TREE: www.iqtree.org/

    """
    
    trees = Phylo.parse(tree, 'newick')
    number = sum(1 for t in trees)
        
    wd = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(msa)))
    args = [exe, '-s', 'msa.fasta', '-zb', '10000', '-au', '-nt', '1', '-pre',
            'autest', '-n', '0']
    if model:
        args.extend(['-m', model])
    else:
        args.extend(['-m', 'TEST'])
    
    shutil.copy(msa, os.path.join(wd, 'msa.fasta'))

    t1, t2 = '', ''
    if number == 1:
        info('One tree found in tree file, inferring ML tree.')
        mltree = mlt(exe, msa, model=model)
        if mltree:
            trees = os.path.join(wd, 'trees.newick')
            if isinstance(tree, str):
                with open(tree) as f1:
                    t1 = f1.readline().rstrip()
            else:
                t1 = tree.format('newick').rstrip()
            with open(trees, 'w') as o, open(mltree) as f2:
                t2 = f2.read().strip()
                o.write('{}\n{}\n'.format(t1, t2))
                args.extend(['-z', 'trees.newick'])
        else:
            error('Infer ML tree failed, AU TEST aborted.')
            sys.exit(1)
    else:
        with open(tree) as f:
            t1, t2 = f.readline().strip(), f.readline().strip()
        shutil.copy(tree, os.path.join(wd, 'trees.newick'))
        args.extend(['-z', 'trees.newick'])
    
    try:
        info('Running AU TEST (IQ-TREE) using the following command:\n\t'
             '{}'.format(' '.join(args)))
        process = Popen(args, cwd=wd, stdout=PIPE, stderr=PIPE,
                        universal_newlines=True)
        code = process.wait()
        if code:
            msg = indent(process.stderr.read(), prefix='\t')
            error('Topology test (AU TEST) failed for {}\n{}'.format(msa, msg))
            sys.exit(1)
        else:
            info('Parsing tree topology test (AU TEST) results.')
            with open(os.path.join(wd, 'autest.iqtree')) as f:
                for line in f:
                    if 'USER TREES' in line:
                        break
        
                for line in f:
                    if line.strip().startswith('1'):
                        aup = float(line.split()[-2])
                        return aup, t1, t2
                error('Parsing tree topology test (AU TEST) failed, no test'
                      'result found.')
                sys.exit(1)
    except OSError:
        error('Invalid IQ-TREE executable (exe) {}, topology test (AU TEST) '
              'failed for {}.'.format(exe, msa))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)


def aut(exe, msa, tree, model='', seed=0, outfile='', verbose=False):
    """
    General use function for performing topology test (AU test).
    
    :param msa: str, path to a FASTA format MSA file.
    :param tree: str, path to a NEWICK format tree file.
    :param exe: str, path to a topology test program's executable.
    :param model: str, name of the model or path of the model file.
    :param seed: int, the seed used to initiate the random number generator.
    :param outfile: str, path to the output file for storing test result.
        If not set, only return the result without save to a file.
    :param verbose: bool, invoke verbose or silent process mode,
        default: False, silent mode.
    :return: tuple, p-value of the AU test (float), first tree (string), and 
        second tree (string).
    """

    logger.setLevel(logging.INFO if verbose else logging.ERROR)

    if os.path.isfile(msa):
        msa = os.path.abspath(msa)
    else:
        error('Alignment {} is not a file or does not exist.'.format(msa))
        sys.exit(1)
    
    if not os.path.isfile(tree):
        error('Tree {} is not a file a does not exist.'.format(tree))
        sys.exit(1)
    
    if not isinstance(model, str):
        error('Argument model accepts a string for the name of the model or '
              'path of the model file')
        sys.exit(1)
        
    try:
        seed = int(seed) if seed else random.randint(0, 1000000)
    except ValueError:
        warn('Invalid seed, generating seed using random number generator.')
        seed = random.randint(0, 1000000)
        
    aup, t1, t2 = _iqtree(exe, msa, tree, model, seed)
    if outfile:
        try:
            with open(outfile, 'w') as o:
                o.write('P-value: {:.4f}\n'
                        'Tree1: {}\n'
                        'Tree2: {}\n'.format(aup, t1, t2))
        except OSError as err:
            error('Save AU TEST result to file failed due to:\n'
                  '\t{}'.format(outfile, err))
    return aup, t1, t2


def main():
    des = 'Perform topology test (AU test) over protein alignments.'
    epilog = """
The minimum requirement for running iMC-aut is an executable of a AU test
program, a multiple sequence alignment (MSA) file (in FASTA format), a tree
file (in NEWICK format).

If you only provide one tree in the tree file, a maximum-likelihood tree will
be inferred and AU-TEST will be performed to test the difference between them.
If a set of trees in NEWICK format was provided in the tree file, only these
trees will be evaluated without reconstructing the ML tree. In both case, only
the p-value of AU TEST for the first tree will be printed out.

Without providing a output file path, AU test results will not be saved to
local file system, but the results will be always printed out.
"""
    formatter = argparse.RawDescriptionHelpFormatter
    parse = argparse.ArgumentParser(description=des, prog='ProtParCon-aut',
                                    usage='%(prog)s EXE MSA TREE [OPTIONS]',
                                    formatter_class=formatter, epilog=epilog)

    parse.add_argument('EXE',
                       help='Path to the executable of the topology test (AU '
                            'test) program.')
    parse.add_argument('MSA',
                       help='Path to the alignment file in FASTA format.')
    parse.add_argument('TREE',
                       help='Path to the tree file or string in NEWICK format.')
    parse.add_argument('-m', '--model',
                       help='Name of the substitution model or filename of the '
                            'model file.')
    parse.add_argument('-o', '--output',
                       help='Path to the output file.')
    parse.add_argument('-s', '--seed',
                       help='The seed for initiating the random number '
                            'generator.')
    parse.add_argument('-v', '--verbose', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')
    
    args = parse.parse_args()
    msa, tree, exe, model = args.MSA, args.TREE, args.EXE, args.model
    
    aup, t1, t2 = aut(exe, msa, tree, model=model, seed=args.seed,
                      outfile=args.out, verbose=args.verbose)
    print(aup)
    print(t1)
    print(t2)


if __name__ == '__main__':
    main()
