#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Common utility functions for various iParCon submodules.
"""

import os
import re
import sys
import shutil
import logging

from io import StringIO
from copy import deepcopy
from collections import namedtuple, defaultdict

from Bio import Phylo, AlignIO

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error

ENDINGS = ['fa', 'fan', 'fas', 'fasta', 'aln', 'clustal', 'clu',
           'phy', 'phylip',  'phylip-relaxed', 'phylip-sequential',
           'stockholm', 'mauve',  'mau', 'emboss', 'emb',
           'nex', 'nexus', 'maf', 'xmfa',
           'newick', 'new',
           'text', 'tsv']
MODEL = namedtuple('model', 'name frequency gamma rates invp type')
AMINO_ACIDS = set('ARNDCQEGHILKMFPSTWYV')


def basename(name):
    """
    Removing file extension, if the extension is in ENDINGS.
    
    :param name: str, a filename.
    :return: str, filename without known extensions for tsv, txt, FASTA, and
        NEWICK as well as some alignment format files.
    """
    
    if '.' in name:
        names = name.split('.')
        if names[-1].lower() in ENDINGS:
            return '.'.join(names[:-1])
        else:
            return name
    else:
        return name
        
        
def modeling(model):
    """
    Parse evolutionary model.

    The model needs to be in the format of MODEL+<FreqType>+<RateType> or a
    model (text) file, where MODEL is a model name (e.g. JTT, WAG, ...),
    FreqType is how the model will handle amino acid frequency (e.g. F, FO,
    or FQ), and RateType is the rate heterogeneity type. Protein mixture
    models are not supported.
    
    :param model: str, name of the model a pathname to the model file.
    :return: namedtuple, in the order of name, frequency, gamma, rates, invp,
        and type.
    """
    
    if os.path.isfile(model):
        return MODEL(os.path.abspath(model), 'empirical', 0, 0, 0, 'custom')
    else:
        name = model.split('+', 1)[0]
        invprop = re.search('\+I[+]*', model)
        invprop = 'estimate' if invprop else 0
        gamma = re.search('\+G([0-9])*\+*', model)
        gamma = int(gamma.group(1)) if gamma else 0
        frequency = re.search('\+(F[OQ]?)\+*', model)
        if frequency:
            fs = {'F': 'empirical', 'FO': 'estimate', 'FQ': 'equal'}
            frequency = fs.get(frequency.group(1), 'empirical')
        else:
            frequency = 'empirical'
        rates = re.search('\+R([0-9])*\+*', model)
        rates = int(rates.group(1)) if rates else 0
        return MODEL(name, frequency, gamma, rates, invprop, 'builtin')
    
    
class Tree(object):
    def __init__(self, tree, leave=False):
        """
        A class handles phylogenetic trees.
        
        :param tree: str, a newick tree file or tree string (must start with
            '(' and end with ';').
        :param leave: bool, whether exit the process or not once an error
            occurred.
        """
        
        if isinstance(tree, str):
            if tree.startswith('(') and tree.endswith(';'):
                tree = StringIO(tree)
            elif os.path.isfile(tree):
                pass
            else:
                error('Invalid tree: {}, tree should be either a NEWICK format '
                      'tree string or tree file.'.format(tree))
                tree = None
        else:
            error('Invalid tree, tree should be a string.')
            tree = None
            
        if tree is None and leave:
            sys.exit(1)
        else:
            tree = Phylo.read(tree, 'newick') if tree else None
            
        if tree:
            leaves = len([c for c in tree.find_clades() if c.is_terminal()])
            nodes = len([c for c in tree.find_clades() if not c.is_terminal()])
            length = sum([c.branch_length if c.branch_length else 0.0
                          for c in tree.find_clades()])
        else:
            leaves, nodes, length = 0, 0, 0.0
            
        self.tree = tree
        self.leaves = leaves
        self.nodes = nodes
        self.length = length
        
    def string(self, ic=True, ic2name=False, nodes=False, brlen=True):
        """
        Get a newick tree string for a tree after manipulating.
        
        :param ic: bool, whether to keep confidence of internal nodes.
        :param ic2name: whether to convert confidence of internal nodes to
            their names.
        :param nodes: bool, discard or keep name of internal nodes.
        :param brlen: bool, discard or keep branch lengths.
        :return: str, a newick tree string.
        """
        
        s = ''
        if self.tree:
            tree = deepcopy(self.tree)
        else:
            return s
        
        n = self.leaves
        for clade in tree.find_clades():
            if ic2name:
                if clade.confidence:
                    if not clade.is_terminal() and not clade.name:
                        clade.name = str(clade.confidence)
                        clade.confidence = None
            if not ic:
                clade.confidence = None
            
            if nodes:
                clade.confidence = None
                if not clade.is_terminal() and not clade.name:
                    clade.name = 'NODE{}'.format(n)
                    n += 1
            else:
                if ic2name:
                    pass
                else:
                    if clade.name and not clade.is_terminal():
                        clade.name = None
                    
            if not brlen:
                clade.branch_length = 0.0
                
        s = tree.format('newick').strip()
        # Remove branch length if brlen set to False
        if not brlen:
            s = re.sub(r':0\.\d+', '', s)
        # Remove branch lengths if only unscaled tree (topology) was given
        if not self.length:
            s = re.sub(r':0\.\d+', '', s)
        return s
        
    def file(self, filename, ic=True, ic2name=False, nodes=False, brlen=True):
        """
        Write a tree to a file after manipulating.
        
        :param filename: str, the name of the tree output file
        :param ic: bool, whether to keep confidence of internal nodes.
        :param ic2name: whether to convert confidence of internal nodes to
            their names.
        :param nodes: bool, discard or keep name of internal nodes.
        :param brlen: bool, discard or keep branch lengths.
        :return: str, a newick tree string.
        """
        s = self.string(ic=ic, ic2name=ic2name, nodes=nodes, brlen=brlen)
        try:
            with open(filename, 'w') as o:
                o.write('{}\n'.format(s))
        except IOError as err:
            error(err)
        return filename


def trim(msa, fmt='fasta', outfile='', verbose=False):
    """
    Remove gaps and ambiguous characters from protein multiple sequence
    alignment (MSA) file or a dictionary object of MSA records.

    :param msa: str or dict, pathname of the protein (MSA) multiple sequence
        alignment file or a dictionary object of of MSA records.
    :param fmt: format of the MSA file, if msa is the pathname of a MSA file.
    :param outfile: pathname for saving trimmed MSA to a file, if not set,
        trimmed sequences will only be returned without saving to a file.
    :param verbose: bool, invoke verbose or silent process mode,
        default: False, silent mode.
    :return:  dict, trimmed MSA in a dictionary.

    .. note::
    
        Ambiguous characters (characters not in ARNDCQEGHILKMFPSTWYV) will be
        treated as gaps and removed. The trimmed alignment file will be saved
        to outfile in FASTA format if outfile is not a empty string and there
        are sites left after trimming. Whether the outfile is set or not,
        `trim()` always returns trimmed MSA in a dictionary (even a empty dict).
        
        Careful with PHYLIP format MSA files, `trim()` use `Bio.AlignIO` to read
        MSA file, users are responsible for given the right name for PHYLIP
        format files, e.g. phylip-relaxed or phylip-sequential.
    """
    
    level = logging.INFO if verbose else logging.ERROR
    logger.setLevel(level)
    
    if isinstance(msa, str):
        if not os.path.isfile(msa):
            error('Alignment {} is not a file or does not exist.'.format(msa))
            sys.exit(1)
        try:
            records = AlignIO.read(msa, fmt)
        except ValueError:
            error('Alignment {} contains non equal length sequences.'.format(
                    msa))
            sys.exit(1)
        length, label = records.get_alignment_length(), msa
        records = {a.id: a.seq for a in records}
    elif isinstance(msa, dict):
        length = list(set([len(v) for v in msa.values()]))
        if not len(length) == 1:
            error('Alignment object contains non equal length sequences.')
            sys.exit(1)
        length, label = length[0], 'alignment record'
        records = msa
    else:
        error('Invalid msa for trimming, msa should be a pathname of MSA file'
              'or a dictionary object of sequence name and sequence pairs.')
        sys.exit(1)
    
    info('Trimming gaps and ambiguous characters for {}'.format(label))
    trimmed = defaultdict(list)
    for i in range(length):
        column = set([s[i] for s in records.values()])
        if column.issubset(AMINO_ACIDS):
            for k, v in records.items():
                trimmed[k].append(v[i])
                
    if all(trimmed.values()):
        trimmed = {k: ''.join(v) for k, v in trimmed.items()}
        if outfile:
            try:
                with open(outfile, 'w') as o:
                    o.writelines('>{}\n{}\n'.format(k, v)
                                 for k, v in trimmed.items())
            except IOError:
                outfile = ''
                warn('IOError, failed to save trimmed alignment.')
    else:
        trimmed = {}
        warn('Successfully removed all gaps and ambiguous characters, but '
             'no site was left after trimming.')
        if outfile:
            warn('No trimmed alignment was saved.')
            outfile = ''
    return trimmed, outfile


def indent(text, prefix, predicate=None):
    """
    A copy of textwrap.indent() method introduce in Python 3.3.
    """
    
    if predicate is None:
        def predicate(line):
            return line.strip()

    def prefixed_lines():
        for line in text.splitlines(True):
            yield (prefix + line if predicate(line) else line)
    return ''.join(prefixed_lines())
        

if __name__ == '__main__':
    pass
