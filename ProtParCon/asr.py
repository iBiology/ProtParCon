#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module is designed to provide a common and easy way to use interface for
ancestral stated reconstruction (ASR) using various ASR programs.

Users are only asked to provide an ASR programs's executable, an aligned 
multiple sequence file (in FASTA format), and a guide tree (in NEWICK) 
format. The general use function ``asr()`` will always return dict object
containing sequence records and a tree object (generated by ``Bio.Phylo`` 
module inside Biopython_) or exit with an error code 1 and an error 
message logged.

Users are recommended only to use function ``asr()`` and avoid to use any
private function inside the module. However, users are strongly recommended
to implement new private functions for additional ASR programs in which they are 
interested and incorporate them into the general use function ``asr()``.

.. _Biopython: https://biopython.org/

"""

import os
import sys
import shutil
import logging
import argparse
import tempfile

from io import StringIO
from subprocess import PIPE, Popen
try:
    from textwrap import indent
except ImportError:
    from ProtParCon.utilities import indent

from Bio import Phylo, AlignIO
from ProtParCon.utilities import basename, modeling, Tree

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error

CODEML_MODELS = os.path.join(os.path.dirname(__file__), 'data')
RAXML_MODELS = ['DAYHOFF', 'DCMUT', 'JTT', 'MTREV', 'WAG', 'RTREV', 'CPREV',
                'VT', 'BLOSUM62', 'MTMAM', 'LG', 'MTART', 'MTZOA', 'PMB',
                'HIVB', 'HIVW', 'JTTD', 'CMUT', 'FLU', 'STMTREV', 'DUMMY',
                'DUMMY2', 'AUTO', 'LG4M', 'LG4X', 'PROT_FILE', 'GTR_UNLINKED',
                'GTR']
FASTML_MODELS = ['JTT', 'LG', 'mtREV', 'cpREV', 'WAG', 'DAYHOFF']

CTL = """   seqfile = {seq}   * sequence data file name
   outfile = mlc  * main result file name
  treefile = {tree}   * tree structure file name

     noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
   verbose = 0   * 1: detailed output, 0: concise output
   runmode = 0   * 0: user tree; 1: semi-automatic; 2: automatic
                 * 3: StepwiseAddition; (4,5): PerturbationNNI; -2: pairwise

   seqtype = 2    * 1: codons; 2: AAs; 3: codons --> AAs
     clock = 0    * 0: no clock, 1: clock, 2: local clock
aaRatefile = {mf}   * only used for aa seqs with model=empirical_(F)

     model = {mn}  * models for AAs or codon-translated AAs:
                    * 0: poisson, 1: proportional, 2: Empirical, 3: Empirical+F
                    * 6: FromCodon, 8: REVaa_0, 9: REVaa(nr=189)

 fix_alpha = 0    * 0: estimate gamma shape parameter; 1: fix it at alpha
     alpha = {alpha}  * initial or fixed alpha, 0: infinity (constant rate)
    Malpha = 0    * different alphas for genes
     ncatG = {gamma}    * # of categories in dG of NSsites models

     getSE = 0    * don't want them, 1: want S.E.s of estimates
RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
"""


def _guess(exe):
    """
    Guess the name of a ancestral states reconstruction (ASR) program according 
    to its executable.

    :param exe: str, path to the executable of an ASR program.
    :return: tuple, name of the ASR program and the corresponding function.
    """
    
    wd = tempfile.mkdtemp()
    try:
        if 'FastML_Wrapper.pl' in exe:
            return 'fastml', _fastml
        process = Popen([exe, '-h'], cwd=wd, stdout=PIPE, stderr=PIPE,
                        universal_newlines=True)
        try:
            outs, errs = process.communicate(timeout=3)
            out = outs or errs
            if 'RAxML' in out:
                return 'raxml', _raxml
            else:
                return 'codeml', _codeml
        except Exception as exc:
            process.terminate()
            process.wait(timeout=10)
            return 'codeml', _codeml
    except OSError:
        error("The exe ({}) is empty or may not be an valid executable of a "
              "ancestral states reconstruction program.".format(exe))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)


def _label(tree, ancestors):
    """
    Relabel internal nodes of a tree and map them to the corresponding name of 
    ancestral sequences.
    
    :param tree: str, a NEWICK format string or file for a tree (must start with
        "(" and end with ';').
    :param ancestors: dict, a dict object stores sequences.
    :return: tuple, a relabeled tree object and a dict object for sequences.
    """
    
    if isinstance(tree, str):
        if os.path.isfile(tree):
            pass
        elif tree.startswith('(') and tree.endswith(';'):
            tree = StringIO(tree)
        else:
            error('Invalid tree encounter, tree relabel aborted.')
            sys.exit(1)
        tree = Phylo.read(tree, 'newick')
    
    number, maps = tree.count_terminals(), {}
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.confidence = None
            number += 1
            old, new = clade.name, 'NODE{}'.format(number)
            maps[old] = new
            clade.name = new
    
    ancestors = {maps.get(k, k): v for k, v in ancestors.items()}
    return tree, ancestors


def _parse(wd):
    """
    Parse the rst file generated by CODEML.
    
    :param wd: str, work directory of  CODEML (inside PAML_ package).
    :return: tuple, a tree object, a dict for sequences, and a list or rates.
    
    .. _PAML: https://www.paml.com/

    """
    
    trees, tree, sequences, ancestors, rates = [], None, [], {}, []
    rst, rs = os.path.join(wd, 'rst'), os.path.join(wd, 'rates')
    if os.path.isfile(rst):
        with open(rst) as f:
            for line in f:
                line = line.strip()
                if line.startswith('(') and line.endswith(';'):
                    trees.append(line)
                elif 'List of extant and reconstructed sequences' in line:
                    break
            
            for line in f:
                line = line.strip()
                if 'Overall accuracy' in line:
                    break
                if line:
                    sequences.append(line)
        
        if trees:
            t1 = Phylo.read(StringIO(trees[0].replace(' ', '')), 'newick')
            brlens = [clade.branch_length for clade in t1.find_clades()]
            
            tree = Phylo.read(StringIO(trees[2].replace(' ', '')), 'newick')
            for i, clade in enumerate(tree.find_clades()):
                clade.branch_length = float(brlens[i]) if brlens[
                    i] else 0.000000
                name = clade.name
                if name and '_' in name:
                    clade.name = name.split('_', maxsplit=1)[1]
                if clade.confidence:
                    clade.name = str(clade.confidence)
                    clade.confidence = None
        else:
            error('Incomplete rst file encounter, no trees were found, parse '
                  'rst file aborted.')
            sys.exit(1)
            
        if sequences:
            for sequence in sequences[1:]:
                blocks = sequence.split()
                if blocks[0] == 'node':
                    name = blocks[1].replace('#', 'NODE')
                    seq = ''.join(blocks[2:])
                else:
                    name, seq = blocks[0], ''.join(blocks[1:])
                ancestors[name] = seq
        else:
            error('Incomplete rst file encounter, no ancestral sequences were '
                  'found, parse rst file aborted.')
            sys.exit(1)
        if tree and ancestors:
            tree, ancestors = _label(tree, ancestors)
            
        if os.path.isfile(rs):
            with open(os.path.join(wd, 'rates')) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        blocks = line.split()
                        if len(blocks) == 5 and blocks[0].isdigit():
                            rates.append(float(blocks[3]))
                        if blocks[0] == 'Site' and rates:
                            break
        else:
            warn('\tParsing PAML result failed (not rates file was found).')
    else:
        error('Parse rst file aborted, the rst file {} does not '
              'exist'.format(rst))
        sys.exit(1)
    return tree, ancestors, rates
    
    
def _write(tree, ancestor, rates, outfile):
    """
    Write tree (object) and ancestor (dict) to a output file.
    
    :param tree: object, tree object.
    :param ancestor: dict, dict object for sequence records.
    :param outfile: str, path to the output file.
    :return: str, path to the output file.
    """
    
    try:
        with open(outfile, 'w') as o:
            o.write('#TREE\t{}\n'.format(tree.format('newick')))
            if rates:
                o.write('#RATES\t{}\n\n'.format('\t'.join([str(r)
                                                         for r in rates])))
            o.writelines('{}\t{}\n'.format(k, v) for k, v in ancestor.items())
    except IOError as err:
        error('Write ancestral reconstruction results to file failed due to:'
              '\n\t{}'.format(err))
        sys.exit(1)
    return outfile


def _codeml(exe, msa, tree, model, gamma, alpha, freq, outfile):
    """
    Reconstruct ancestral sequences using CODEML (inside PAML_ package).
    
    :param exe: str, path to the executable of an ASR program.
    :param msa: str, path to the MSA file (must in FASTA format).
    :param tree: str, path to the tree file (must in NEWICK format) or a NEWICK
        format tree string (must start with "(" and end with ";").
    :param model: namedtuple, substitution model for ASR.
    :param gamma: int, The number of categories for the discrete gamma rate
        heterogeneity model. 
    :param freq: str, the equilibrium frequencies of the twenty amino acids.
    :param alpha: float, the shape (alpha) for the gamma rate heterogeneity.
    :param outfile: str, path to the output file. 

    :return: tuple, a tree object, a dict for sequences, and a list or rates.
    
    .. note::
        See doc string of function ``asr()`` for details of all arguments.
    
    .. _PAML: https://www.paml.com/
    
    """
    
    cwd = os.getcwd()
    if model.type == 'custom':
        mf = model.name
        info('Use custom model file {} for ancestral states '
             'reconstruction.'.format(mf))
    else:
        name = model.name
        info('Use {} model for ancestral states reconstruction.'.format(name))
        mf = os.path.join(CODEML_MODELS, '{}'.format(name.lower()))
        if not os.path.isfile(mf):
            error('Failed to find model file for model {}, aborted.'.format(
                    name))
            sys.exit(1)
    
    wd, tf = tempfile.mkdtemp(dir=os.path.dirname(msa)), 'codeml.tree.newick'
    tf = tree.file(os.path.join(wd, tf), brlen=False)
        
    parameters = {'seq': msa, 'tree': tf, 'mf': mf}
    if model.frequency == 'estimate':
        parameters['mn'] = 3
    else:
        parameters['mn'] = 2
    parameters['alpha'] = alpha if alpha else 0.5
    gamma = gamma or model.gamma
    parameters['gamma'] = gamma if gamma else 4
    
    with open(os.path.join(wd, 'ctl.dat'), 'w') as o:
        o.write(CTL.format(**parameters))

    try:
        # No clue why stdout and stderr PIPEs are here (not in _raxml()) will cause
        # the following warn in unittest:
        #
        # ResourceWarning: unclosed file <_io.TextIOWrapper name=3
        # encoding='cp1252'>
        #
        # ResourceWarning: unclosed file <_io.TextIOWrapper name=4
        # encoding='cp1252'>
        info('Reconstructing ancestral states for {} using CODEML.'.format(msa))
        process = Popen([exe, 'ctl.dat'], cwd=wd, stdout=PIPE,
                        stderr=PIPE, universal_newlines=True)

        code = process.wait()
        msg = process.stdout.read() or process.stderr.read()
        if code:
            error('Ancestral reconstruction via CODEML failed for {} due to:'
                  '\n{}'.format(msa, indent(msg, prefix='\t')))
            sys.exit(1)
        else:
            info('Parsing ancestral sequence reconstruction results.')
            tree, ancestors, rates = _parse(wd)
            
            outfile = _write(tree, ancestors, rates, outfile)
            info('Successfully save ancestral states reconstruction '
                 'results to {}.'.format(outfile))
            return outfile
    except OSError:
        error('Invalid PAML (CODEML) executable {}, running CODEML failed for '
              '{}.'.format(exe, msa))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)
        
        
def _raxml(exe, msa, tree, model, gamma, alpha, freq, outfile):
    """
    Reconstruct ancestral sequences using RAxML_.
    
    :param exe: str, path to the executable of an ASR program.
    :param msa: str, path to the MSA file (must in FASTA format).
    :param tree: str, path to the tree file (must in NEWICK format) or a NEWICK
        format tree string (must start with "(" and end with ";").
    :param model: namedtuple, substitution model for ASR.
    :param gamma: int, The number of categories for the discrete gamma rate
        heterogeneity model.
    :param freq: str, the equilibrium frequencies of the twenty amino acids.
    :param alpha: float, the shape (alpha) for the gamma rate heterogeneity.
    :param outfile: str, path to the output file.

    :return: tuple, a tree object and a dict for sequences.
    
    .. note::
        See doc string of function ``asr()`` for details of all arguments.
        
    .. _RAxML: https://sco.h-its.org/exelixis/web/software/raxml/
    """
    
    if model.type == 'custom':
        mf = model.name
        name = 'WAG'
        info('Use model file {} for ancestral states reconstruction.'.format(mf))
    else:
        name = model.name
        if name.upper() in RAXML_MODELS:
            mf = ''
            info('Use {} model for ancestral states '
                 'reconstruction.'.format(name))
        else:
            error('RAxML does not accept {} model, aborted.'.format(name))
            sys.exit(1)
            
    wd, tf = tempfile.mkdtemp(dir=os.path.dirname(msa)), 'raxml.tree.newick'
    tf = tree.file(os.path.join(wd, tf), brlen=False)
        
    cmd = [exe, '-f', 'A', '-s', msa, '-t', tf, '-n', 'iMC']
    m = 'PROTGAMMA' if (gamma or model.gamma) else 'PROTCAT'
    m += name.upper()
    freq = 'F' if freq == 'estimate' or model.frequency == 'estimate' else 'X'
    if 'AUTO' not in m:
        m += freq
    cmd.extend(['-m', m])
    if mf:
        cmd.extend(['-P', mf])
    cmd.append('--silent')
    
    try:
        info('Reconstructing ancestral states for {} using RAxML.'.format(msa))
        process = Popen(cmd, cwd=wd, stdout=PIPE,
                        stderr=PIPE, universal_newlines=True)
        code = process.wait()
        msg = process.stdout.read() or process.stderr.read()
        # Sometime RAxML does not return a non-zero code when it runs error
        if code:
            error('Ancestral reconstruction via RAxML failed for {} due to:'
                  '\n{}'.format(msa, indent(msg, prefix='\t')))
            sys.exit(1)
        else:
            ancestor = os.path.join(wd, 'RAxML_marginalAncestralStates.iMC')
            # Check again to see if reconstruction success
            if not os.path.isfile(ancestor):
                msg = '\n'.join([line for line in msg.splitlines()
                                 if line.strip().startswith('ERROR')])
                error('Ancestral reconstruction via RAxML failed for {} due to:'
                      '\n{}'.format(msa, indent(msg, prefix='\t')))
                sys.exit(1)
            info('Parsing ancestral sequence reconstruction results.')
            with open(ancestor) as handle:
                ancestor = dict(line.strip().split() for line in handle)
            tree = os.path.join(wd, 'RAxML_nodeLabelledRootedTree.iMC')
            tree = Phylo.read(tree, 'newick')
            for clade in tree.find_clades():
                if clade.confidence and not clade.name:
                    clade.name = str(clade.confidence)
            tree, ancestor = _label(tree, ancestor)
            for record in AlignIO.read(msa, 'fasta'):
                ancestor[record.id] = record.seq
            
            _write(tree, ancestor, [], outfile)
            info('Successfully save ancestral states reconstruction '
                 'results to {}.'.format(outfile))
            return outfile
    except OSError as err:
        print(err)
        error('Invalid RAxML executable {}, running RAxML failed for '
              '{}.'.format(exe, msa))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)
        
        
def _fastml(exe, msa, tree, model, gamma, alpha, freq, outfile):
    """
    Reconstruct ancestral sequences using FastML_.

    :param exe: str, path to the executable of an ASR program.
    :param msa: str, path to the MSA file (must in FASTA format).
    :param tree: str, path to the tree file (must in NEWICK format) or a NEWICK
        format tree string (must start with "(" and end with ";").
    :param model: namedtuple, substitution model for ASR.
    :param gamma: int, The number of categories for the discrete gamma rate
        heterogeneity model.
    :param freq: str, the equilibrium frequencies of the twenty amino acids.
    :param alpha: float, the shape (alpha) for the gamma rate heterogeneity.
    :param outfile: str, path to the output file.

    :return: tuple, a tree object and a dict for sequences.

    .. note::
        See doc string of function ``asr()`` for details of all arguments.

    .. _FastML: http://fastml.tau.ac.il/
    """
    
    if model.type == 'custom':
        mf = model.name
        name = 'JTT'
        info('Use model file {} for ancestral states reconstruction.'.format(mf))
    else:
        name = model.name
        if name in FASTML_MODELS:
            mf = ''
            info('Use {} model for ancestral states reconstruction.'.format(name))
        else:
            error('RAxML does not accept {} model, aborted.'.format(name))
            sys.exit(1)
    
    wd, tf = tempfile.mkdtemp(dir=os.path.dirname(msa)), 'fastml.tree.newick'
    tf = tree.file(os.path.join(wd, tf), brlen=False)
    
    cmd = ['perl', exe, '--MSA_File', msa, '-seqType', 'AA', '--Tree', tf,
           '--outDir', wd]
    
    if gamma or model.gamma:
        cmd.extend(['--UseGamma', 'yes'])
        if alpha:
            cmd.extend(['--Alpha', str(alpha)])
    
    try:
        info('Reconstructing ancestral states for {} using FastML.'.format(msa))
        process = Popen(cmd, cwd=wd, stdout=PIPE, stderr=PIPE,
                        universal_newlines=True)
        code = process.wait()
        msg = process.stdout.read() or process.stderr.read()
        if code:
            error('Ancestral reconstruction via FastML failed for {} due to:'
                  '\n{}'.format(msa, indent(msg, prefix='\t')))
            sys.exit(1)
        else:
            ancestor = os.path.join(wd, 'seq.marginal.txt')
            # Check again to see if reconstruction success
            if not os.path.isfile(ancestor):
                msg = '\n'.join([line for line in msg.splitlines() if
                                 line.strip().startswith('ERROR')])
                error('Ancestral reconstruction via RAxML failed for {} due to:'
                      '\n{}'.format(msa, indent(msg, prefix='\t')))
                sys.exit(1)
            info('Parsing ancestral sequence reconstruction results.')
            ancestor = {a.id: a.seq for a in AlignIO.read(ancestor, 'fasta')}

            tree = os.path.join(wd, 'tree.newick.txt')
            tree = Phylo.read(tree, 'newick')

            _write(tree, ancestor, [], outfile)
            info('Successfully save ancestral states reconstruction '
                 'results to {}.'.format(outfile))
            return outfile
    except OSError as err:
        error('Invalid FastML executable {}, running RAxML failed for '
              '{}.'.format(exe, msa))
        sys.exit(1)
    finally:
        shutil.rmtree(wd)
    
        
def asr(exe, msa, tree, model, gamma=4, alpha=1.8, freq='',
        outfile='', verbose=False):
    """
    General use function for (marginal) ancestral states reconstruction (ASR).

    :param exe: str, path to the executable of an ASR program.
    :param msa: str, path to the MSA file (must in FASTA format).
    :param tree: str, path to the tree file (must in NEWICK format) or a NEWICK
        format tree string (must start with "(" and end with ";").
    :param model: str, substitution model for ASR. Either a path to a model
        file or a valid model string (name of an empirical model plus some 
        other options like gamma category and equilibrium frequency option).
        If a model file is in use, the file format of the model file depends
        on the ASR program, see the its documentation for details.
    :param gamma: int, The number of categories for the discrete gamma rate
        heterogeneity model. Without setting gamma, RAxML will use CAT model
        instead, while CODEML will use 4 gamma categories.
    :param freq: str, the base frequencies of the twenty amino acids.
        Accept empirical, or estimate, where empirical will set frequencies
        use the empirical values associated with the specified substitution
        model, and estimate will use a ML estimate of base frequencies.
    :param alpha: float, the shape (alpha) for the gamma rate heterogeneity.
    :param outfile: str, path to the output file. Whiteout setting, results
        of ancestral states reconstruction will be saved using the filename
        `[basename].[asrer].tsv`, where basename is the filename of MSA file
        without known FASTA file extension, asrer is the name of the ASR
        program (in lower case). The first line of the file will start with
        '#TREE' and followed by a TAB (\t) and then a NEWICK formatted tree
        string, the internal nodes were labeled. The second line of the tsv
        file is intentionally left as a blank line and the rest lines of the
        file are tab separated sequence IDs and amino acid sequences.
    :param verbose: bool, invoke verbose or silent (default) process mode.
    :return: tuple, the paths of the ancestral states file.

    .. note::
    
        If a tree (with branch lengths and/or internal nodes labeled) is
        provided, the branch lengths and internal node labels) will be ignored.
        
        If the model name combined with Gamma category numbers, i.e. JTT+G4,
        WAG+G8, etc., only the name of the model will be used. For all models
        contain G letter, a discrete Gamma model will be used to account for
        among-site rate variation. If there is a number after letter G,
        the number will be used to define number of categories in CODEML. For
        RAxML, the number of categories will always be set to 4 if G presented.
        
    """
    
    level = logging.INFO if verbose else logging.ERROR
    logger.setLevel(level)
    
    if os.path.isfile(msa):
        msa = os.path.abspath(msa)
    else:
        error('Ancestral reconstruction aborted, msa {} is not a file or '
              'does not exist.'.format(msa))
        sys.exit(1)
    
    tree = Tree(tree, leave=True)
        
    if not isinstance(model, str):
        error('Ancestral reconstruction aborted, model {} is not a valid '
              'model name or model file.'.format(model))
        sys.exit(1)
        
    model = modeling(model)
    asrer, func = _guess(exe)
    if not outfile or outfile == 'iMC-default':
        outfile = '{}.{}.tsv'.format(basename(msa), asrer)

    out = func(exe, msa, tree, model, gamma, alpha, freq, outfile)
    return out


def main():
    des = 'Ancestral states reconstruction (ML method) over protein alignments.'
    epilog = """
The minimum requirement for running iMC-aut is an executable of a ancestral
states reconstruction program, a multiple sequence alignment (MSA) file in FASTA
format, and a guide tree (or topology) file in NEWICK format file or string. 

If a tree (with branch lengths and/or internal nodes labeled) was provided, the
branch lengths and internal node labels will be ignored. All branch length will
be estimated during states reconstruction and all names of internal nodes will
be correspondingly relabeled to match the sequence IDs.

If the model name combined with Gamma category numbers, i.e. JTT+G4, WAG+G8, ...
a discrete Gamma model will be used to account for among-site rate variation. If
there is a number after letter G, the number will be used to define number of
categories in CODEML. Without the number, argument gamma will be checked, if not
set, 4 categories will be applied to CODEML and RAxML. For RAxML, if no G in the
model and no gamma was set, CAT model will be used.

Ancestral states will be always saved to a tab separated text file. In case
of no output filename was set via -o flag, a filename `[basename].[asrer].tsv`
will be used, where basename is the filename of MSA file without known FASTA
file extension, asrer is the name of the ASR program (in lower case). The first
line of the file will start with  '#TREE' and followed by a TAB ('\t') and then
a NEWICK formatted tree string, the internal nodes were labeled. The second line
of the tsv file is intentionally left as a blank line and the rest lines of the
file are tab separated sequence IDs and amino acid sequences.
"""
    
    formatter = argparse.RawDescriptionHelpFormatter
    parse = argparse.ArgumentParser(description=des, prog='ProtParCon-anc',
                                    usage='%(prog)s EXE MSA TREE [OPTIONS]',
                                    formatter_class=formatter, epilog=epilog)
    
    parse.add_argument('EXE',
                       help='Pathname of the executable of the ancestral '
                            'states reconstruction program.')
    parse.add_argument('MSA',
                       help='Pathname of the multiple sequence alignment (MSA) '
                            'file in FASTA format.')
    parse.add_argument('TREE',
                       help='Pathname of the guide tree (or topology) file '
                            'in NEWICK format.')
    parse.add_argument('-model', default='JTT',
                       help='Name of the evolutionary model or filename of the '
                            'model file.')
    parse.add_argument('-gamma',
                       help='The number of categories for the discrete gamma '
                            'rate heterogeneity model.')
    parse.add_argument('-alpha',
                       help='The shape (alpha) for the gamma rate '
                            'heterogeneity.')
    parse.add_argument('-freq',
                       help='Comma separated state frequencies.')
    parse.add_argument('-o',
                       help='Path of the output file.')
    parse.add_argument('-v', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')
    
    args = parse.parse_args()
    exe, tree, msa, model = args.EXE, args.TREE, args.MSA, args.model
    freq, gamma, alpha, v = args.freq, args.gamma, args.alpha, args.v
    out = args.o

    outfile = out if out else 'iMC-default'
    
    asr(exe, msa, tree, model, gamma=gamma, alpha=alpha, freq=freq,
        verbose=v, outfile=outfile)


if __name__ == '__main__':
    main()
