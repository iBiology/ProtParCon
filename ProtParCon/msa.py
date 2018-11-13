#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module is designed to provide a common and easy way to use interface for
aligning multiple sequences using various align programs.

Users are only asked to provide an alignment programs's executable and an 
multiple sequence file in FASTA format. The general use function ``msa()`` 
will always return the pathname of the alignment output file or exit with 
an error code 1 and an error message logged.

Users are recommended to only use function ``msa()`` and avoid to use any
private functions inside the module. However, we strongly recommend users 
to implement new private functions for additional alignment programs in which 
they are interested and incorporate them into the general use function 
``msa()``.
"""

import os
import sys
import shutil
import logging
import argparse
import tempfile
try:
    from textwrap import indent
except ImportError:
    from ProtParCon.utilities import indent
from subprocess import PIPE, Popen

from ProtParCon.utilities import basename, trim

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error


def _guess(exe):
    """
    Guess the name of a multiple sequence alignment (MSA) program according to 
    its executable.
    
    :param exe: str, path to the executable of a MSA program.
    :return: tuple, name of the MSA program and the corresponding function.
    """
    
    aligner, func = None, None
    aligners = {'muscle': _muscle, 'mafft': _mafft,
                'clustal': _clustal, 't-coffee': _tcoffee}
    try:
        process = Popen([exe, '--help'], stdout=PIPE, stderr=PIPE,
                        universal_newlines=True)
        outs, errs = process.communicate(timeout=10)
        out = outs or errs
        out = out[:100]
        for a, f in aligners.items():
            if a.upper() in out or a.title() in out:
                aligner, func = a, f
                break
    except OSError:
        error('The aligner exe: {} is empty or may not be an valid executable '
              'of a sequence align program.'.format(exe))
        sys.exit(1)
    return aligner, func


def _mafft(exe, seq, outfile):
    """
    Align multiple sequences using MAFFT_.
    
    :param exe: str, path to the executable of a multiple sequence align
        program.
    :param seq: str, path to the multiple sequence file (must in FASTA format).
    :param seq: str, path to the aligned sequence output file (in FASTA format).
    
    :return: str, path to the aligned sequence output file (in FASTA format).
    
    .. _MAFFT: https://mafft.cbrc.jp/alignment/software/
    """
    
    args = [exe, '--quiet', seq]
    try:
        with open(outfile, 'w') as stdout:
            process = Popen(args, stdout=stdout, stderr=PIPE,
                            universal_newlines=True)
            code = process.wait()
    except OSError:
        msg = 'Failed to write alignment to outfile {}'.format(outfile)
        error('Aligning sequence: {} via MUSCLE failed due to:\n\tIOError, '
              'failed to write alignment to outfile: {}.'.format(seq, outfile))
        sys.exit(1)
    if code:
        if os.path.isfile(outfile):
            os.remove(outfile)
        msg = indent(process.stderr.read(), prefix='\t')
        process.stderr.close()
        error('Aligning sequence: {} via MAFFT failed due to:\n{}.'.format(seq,
              msg))
        sys.exit(1)
    
    return outfile


def _muscle(exe, seq, outfile):
    """
    Align multiple sequences using MUSCLE_.
    
    :param exe: str, path to the executable of a multiple sequence align
        program.
    :param seq: str, path to the multiple sequence file (must in FASTA format).
    :param seq: str, path to the aligned sequence output file (in FASTA format).
    
    :return: str, path to the aligned sequence output file (in FASTA format).
    
    .. _MUSCLE: https://www.drive5.com/muscle/
    """
    
    args = [exe, '-in', seq, '-out', outfile, '-quiet']
    process = Popen(args, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    code = process.wait()
    if code:
        if os.path.isfile(outfile):
            os.remove(outfile)
        msg = process.stderr.read() or process.stdout.read()
        process.stdout.close(), process.stderr.close()
        error('Aligning sequence: {} via MUSCLE failed due to:\n{}.'.format(
              seq, indent(msg, prefix='\t')))
        sys.exit(1)
    else:
        process.poll()
    
    return outfile


def _clustal(exe, seq, outfile):
    """
    Align multiple sequences using `CLUSTAL (OMEGA)`_.
    
    :param exe: str, path to the executable of a multiple sequence align
        program.
    :param seq: str, path to the multiple sequence file (must in FASTA format).
    :param seq: str, path to the aligned sequence output file (in FASTA format).
    
    :return: str, path to the aligned sequence output file (in FASTA format).
    
    .. _`CLUSTAL (OMEGA)`: http://www.clustal.org/omega/
    """
    
    args = [exe, '-i', seq, '-o', outfile]
    if os.path.isfile(outfile):
        args.append('--force')
    process = Popen(args, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    code = process.wait()
    if code:
        if os.path.isfile(outfile):
            os.remove(outfile)
        msg = process.stderr.read() or process.stdout.read()
        process.stdout.close(), process.stderr.close()
        error('Aligning sequence: {} via CLUSTAL failed due to:\n{}.'.format(
              seq, indent(msg, prefix='\t')))
        sys.exit(1)
    
    return outfile


def _tcoffee(exe, seq, outfile):
    """
    Align multiple sequences using `T-COFFEE`_.

    :param exe: str, path to the executable of a multiple sequence align
        program.
    :param seq: str, path to the multiple sequence file (must in FASTA format).
    :param seq: str, path to the aligned sequence output file (in FASTA format).

    :return: str, path to the aligned sequence output file (in FASTA format).

    .. _`T-COFFEE`: http://www.tcoffee.org/Projects/tcoffee/
    """
    
    wd = tempfile.mkdtemp(dir=os.path.dirname(seq))
    args = [exe, '-in', seq, '-outfile', outfile, '-output', 'fasta_aln',
            '-run_name', 't-coffee-alignment', '-quiet']
    try:
        process = Popen(args, stderr=PIPE, universal_newlines=True, cwd=wd)
        code = process.wait()
        if code:
            if os.path.isfile(outfile):
                os.remove(outfile)
            msg = process.stderr.read() or process.stdout.read()
            process.stderr.close()
            error('Aligning sequence: {} via T-COFFEE failed due to:\n{}.'
                  .format(seq, indent(msg, prefix='\t')))
            sys.exit(1)
    finally:
        shutil.rmtree(wd)
    
    return outfile


def msa(exe, seq, outfile='', verbose=False, trimming=False):
    """
    General use function for multiple sequence alignment (MSA).

    :param exe: str, path to the executable of a MSA program.
    :param seq: str, path to the multiple sequence file (must in FASTA format).
    :param outfile: str, path to the aligned sequence output (FASTA) file,
        default: [basename].[aligner].fasta, where basename is the filename of 
        the sequence file without known FASTA file extension, aligner is the 
        name of the aligner program in lowercase, and fasta is the extension 
        for fasta format file.
    :param verbose: bool, invoke verbose or silent process mode,
        default: False, silent mode.
    :param trimming: bool, trim gaps and ambiguous sites if True, otherwise,
        leave them untouched.
    :return: str, path to the aligned sequence output file (in FASTA format).
    """
    
    level = logging.INFO if verbose else logging.ERROR
    logger.setLevel(level)
    
    if os.path.isfile(seq):
        sequence = os.path.abspath(seq)
        
        if exe:
            aligner, func = _guess(exe)
            if func is None:
                error('Invalid or unsupported aligner executable (exe): {}, '
                      'alignment aborted.'.format(exe))
                sys.exit(1)
        else:
            error('Invalid aligner executable (exe), empty string, sequence '
                  'alignment aborted.')
            sys.exit(1)
        
        if not outfile:
            outfile = '.'.join([basename(sequence), aligner, 'fasta'])

        info('Aligning sequence {} using {}.'.format(sequence, aligner.upper()))
        outfile = func(exe, sequence, outfile)
        info('Successfully align sequence, alignment was saved to '
             '{}.'.format(outfile))
    else:
        error('Sequence: {} is not a file or does not exist.'.format(seq))
        sys.exit(1)
    if trimming:
        clean = ''.join([basename(outfile), '.trimmed.fa'])
        if os.path.isfile(clean):
            info('Found pre-existing trimmed alignment.')
        else:
            trim(outfile, outfile=clean, verbose=verbose)
            outfile = clean
    return outfile


def main():
    des = 'Common interface for aligning multiple sequences.'
    epilog = """
In order to make this interface common and simple, only FASTA format is the
acceptable sequence file format.

Without specifying the alignment output, alignment will be saved to a file
named in the format of [basename].[aligner].fasta, where basename is the
filename of the sequence file without extension, aligner is the name of the
align program (in lower case), and fasta is the extension for FASTA format 
file.

Under silent process mode, only errors was logged, while under verbose mode,
all errors, warnings and information about processing details will be logged.
"""

    formatter = argparse.RawDescriptionHelpFormatter
    parse = argparse.ArgumentParser(description=des, prog='ProtParCon-msa',
                                    usage='%(prog)s EXECUTABLE SEQUENCE',
                                    formatter_class=formatter,
                                    epilog=epilog)
    
    parse.add_argument('EXECUTABLE',
                       help='path to the executable of a multiple sequence '
                            'align program.')
    parse.add_argument('SEQUENCE',
                       help='Path to the multiple sequence input file '
                            '(must in FASTA format).')
    parse.add_argument('-o',
                       help='Path to the aligned multiple sequence output'
                            'file (in FASTA format).')
    parse.add_argument('-t', action='store_true',
                       help='Trim gaps and ambiguous sites or leave them '
                            'untouched (default).')
    parse.add_argument('-v', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')
    
    args = parse.parse_args()
    seq, exe, out, verbose = args.SEQUENCE, args.EXECUTABLE, args.o, args.v
    trimming = args.t
    msa(exe, seq, outfile=out, verbose=verbose, trimming=trimming)


if __name__ == '__main__':
    pass
