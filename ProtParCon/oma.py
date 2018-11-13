#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import binascii
import argparse
import collections

from itertools import chain
from urllib.request import urlretrieve, urlcleanup

from Bio import SeqIO, SeqRecord

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iMC]')
warn, info, error = logger.warning, logger.info, logger.error

OMA_SPECIES = 'https://omabrowser.org/All/oma-species.txt'
OMA_GROUPS = 'https://omabrowser.org/All/oma-groups.txt.gz'
OMA_SEQUENCES = 'https://omabrowser.org/All/oma-seqs.fa.gz'


def _download(url):
    """
    Private function for checking and downloading data files from OMA database
    if necessary.
    """
    
    filename = url.split('/')[-1]
    if os.path.isfile(filename):
        info('Using pre-existed file {} from local system.'.format(filename))
    else:
        info('Downloading {} from OMA Database.'.format(url.split('/')[-1]))
        filename, _ = urlretrieve(url, filename)
    return filename


def _gz(filename):
    """
    Private function for checking if the file is gzipped file.
    """
    
    with open(filename, 'rb') as f:
        return binascii.hexlify(f.read(2)) == b'1f8b'


def _lines(filename):
    """
    Private function for parsing tab separated lines of regular/gzip text
    file into blocks.
    """
    
    handle = gzip.open(filename, 'rt') if _gz(filename) else open(filename)
    for line in handle:
        if not line.startswith('#'):
            yield line.strip().split('\t')


def _group(codes, group_file):
    """
    Private function for retrieving orthologous groups.
    """
    
    groups, size = {}, len(codes)
    group_temp = 'oma_temporary_groups.tsv'
    if os.path.isfile(group_temp):
        info('Loading pre-existed temporary OMA ortholog groups (oma_temporary_'
             'groups.tsv) ...')
        for blocks in _lines(group_temp):
            groups[blocks[0]] = blocks[1:]
    else:
        info('Parsing OMA ortholog groups (oma-groups.txt.gz) ...')
        for blocks in _lines(group_file):
            number, finger, entries = blocks[0], blocks[1], blocks[2:]
            ids = [entry for entry in entries if entry[:5] in codes]
            if size == len(set(i[:5] for i in ids)):
                groups[finger] = ids
        if groups:
            with open(group_temp, 'w') as o:
                o.writelines('{}\t{}\n'.format(k, '\t'.join(v))
                             for k, v in groups.items())
    info('Yield {} one-to-one ortholog groups for {} query items.'.format(
        len(groups), size))
    return groups


def _seq(codes, seq_file):
    """
    Private function for retrieving orthologous protein sequences.
    """
    
    seq_temp = 'oma_temporary_sequences.fasta'
    if os.path.isfile(seq_temp):
        info('Indexing pre-existed temporary protein sequences ('
             'oma_temporary_sequences.fasta) ... ')
        seqs = SeqIO.index(seq_temp, 'fasta')
    else:
        info('Parsing OMA protein sequences (oma-seqs.fa.gz) ... ')
        handle = gzip.open(seq_file, 'rt') if _gz(seq_file) else open(seq_file)
        records = SeqIO.parse(handle, 'fasta')
        seqs = {record.id: record for record in records if
                record.id[:5] in codes}
        SeqIO.write(seqs.values(), seq_temp, 'fasta')
        handle.close()
    return seqs


def oma(query, group='', seq='', outdir='', verbose=False):
    """
    Query and retrieve 1:1 orthologous protein sequences from OMA
    orthology Database (https://omabrowser.org)

    :param query: sequence object consists of taxa IDs, OMA codes or
        scientific names or mix of them.
        
        For example (all these four queries will yield the same result):
            [9606,9913,9823,9685],
            
            ('Homo sapiens', 'Bos taurus',  'Sus scrofa', 'Felis catus'),
            
            ['HUMAN', 'BOVIN', 'PIGXX', 'FELCA'],
            
            [9606, 'BOVIN', 'Sus_scrofa', 'Felis_catus']
    :param group: str, path to the OMA groups gzip (or unzipped) file
        stored on local drive.
    :param seq: str, path to the OMA sequences gzip (or unzipped) file
        stored on local drive.
    :param outdir: str, output directory, if not set, current work directory
        will be used.
    :param verbose: bool, invoke verbose or silent process mode, default:
        False, silent mode.
    :return: list, consists of absolute paths of protein sequences saved for
        all one-to-one orthologous groups.

    .. notes::
        Without assigning group and/or seq, the corresponding file will be
        automatically downloaded from OMA Database.
        
        Orthologous protein sequence file will be saved in FASTA format using
        [finger].fasta as the file name.
    """
    
    level = logging.INFO if verbose else logging.ERROR
    logger.setLevel(level)
    
    if isinstance(query, (list, tuple, collections.Iterable)):
        queries = []
        for q in query:
            try:
                queries.append(str(q).replace('_', ' '))
            except ValueError:
                error('Invalid query item found: {}.'.format(q))
                error('Query items accept taxa ids (integer or string) or '
                      'scientific names (string) or mix of them.')
                sys.exit(1)
    else:
        error('Invalid query, query only accepts Python sequence object.')
        sys.exit(1)
    
    if os.path.isdir(outdir):
        cwd = os.path.abspath(outdir)
    else:
        try:
            os.mkdir(outdir)
            cwd = os.path.abspath(outdir)
        except IOError:
            cwd = os.getcwd()
    os.chdir(cwd)
    
    size = len(queries)
    info('Searching OMA database with {} query items ...'.format(size))
    
    species = 'oma-species.txt' if os.path.isfile(
        'oma-species.txt') else _download(OMA_SPECIES)
    codes = {}
    for blocks in _lines(species):
        code, taxid, name = blocks[:3]
        if (name in queries) or (taxid in queries) or (code in queries):
            codes[code] = name.replace(' ', '_')
    
    if len(codes) != size:
        records = list(chain.from_iterable(
                [(k, v.replace('_', ' ')) for k, v in codes.items()]))
        absent = [q for q in query if q not in records]
        error('The following {} query items are not found in OMA:'
              '\n\t{}.'.format(len(absent), ', '.join(absent)))
        sys.exit(1)
    
    group_file = group if os.path.isfile(group) else _download(OMA_GROUPS)
    groups = _group(codes, group_file)
    
    sequences = []
    if groups:
        size = len(groups)
        seq_file = seq if os.path.isfile(seq) else _download(OMA_SEQUENCES)
        
        seqs = _seq(codes, seq_file)
        for i, (finger, entries) in enumerate(groups.items(), start=1):
            filename = os.path.join(cwd, '{}.fasta'.format(finger))
            if os.path.isfile(filename):
                sequences.append(filename)
            else:
                sequence = [seqs.get(entry, None) for entry in entries]
                if all(sequence):
                    sequence = [
                        SeqRecord.SeqRecord(s.seq, codes[s.id[:5]], '', '') for
                        s in sequence]
                    SeqIO.write(sequence, filename, 'fasta')
                    sequences.append(filename)
                else:
                    warn('Missing sequence was found in ortholog group: {}, '
                         'skipped.'.format(finger))
        info('Yield {} orthologous protein sequences for {} query '
             'items.'.format(len(sequences), size))
    return sequences


def main():
    des = 'Query and retrieve orthologous proteins from OMA orthology Database.'
    epilog = """
If the query string contains scientific names, space needs to be
replaced with underscore ('_').

Without assigning group and/or seq, the corresponding file will be
automatically downloaded from OMA Database.

Orthologous protein sequence file will be saved in FASTA format using
[finger].fasta as the file name.
    """
    
    formatter = argparse.RawDescriptionHelpFormatter
    parse = argparse.ArgumentParser(description=des, prog='ProtParCon-oma',
                                    usage='%(prog)s QUERY [OPTIONS]',
                                    formatter_class=formatter, epilog=epilog)
    
    parse.add_argument('QUERY',
                       help='Comma separated string consists of taxa IDs, '
                            'OMA codes or scientific names or mix of them.')
    parse.add_argument('-g',
                       help='Path to the OMA groups gzip (or unzipped) file.')
    parse.add_argument('-s',
                       help='Path to the OMA sequences gzip (or unzipped) '
                            'file.')
    parse.add_argument('-o', help='Path of the output directory.')
    parse.add_argument('-v', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')
    
    args = parse.parse_args()
    query, g, s, o, v = args.QUERY, args.g, args.s, args.o, args.v
    query = query.split(',')
    oma(query, group=g, seq=s, outdir=o, verbose=v)


if __name__ == '__main__':
    pass
    # folder = r'C:\Users\tianz\Downloads\iparcon-dataset\oma'
    # taxa = [('ORNAN', 9258, 'Ornithorhynchus anatinus'),
    #         ('TURTR', 9739, 'Tursiops truncatus'),
    #         ('SARHA', 9305, 'Sarcophilus harrisii'),
    #         ('FELCA', 9685, 'Felis catus'),
    #         ('MYOLU', 59463, 'Myotis lucifugus'),
    #         ('AILME', 9646, 'Ailuropoda melanoleuca'),
    #         ('CANLF', 9615, 'Canis lupus familiaris'),
    #         ('OTOGA', 30611, 'Otolemur garnettii'),
    #         ('HORSE', 9796, 'Equus caballus'),
    #         ('LOXAF', 9785, 'Loxodonta africana'),
    #         ('BOVIN', 9913, 'Bos taurus'), ('SHEEP', 9940, 'Ovis aries'),
    #         ('PIGXX', 9823, 'Sus scrofa'),
    #         ('GORGO', 9595, 'Gorilla gorilla gorilla'),
    #         ('CALJA', 9483, 'Callithrix jacchus'),
    #         ('HUMAN', 9606, 'Homo sapiens'),
    #         ('DASNO', 9361, 'Dasypus novemcinctus'),
    #         ('MACMU', 9544, 'Macaca mulatta'),
    #         ('RATNO', 10116, 'Rattus norvegicus'),
    #         ('MOUSE', 10090, 'Mus musculus')]
    #
    # t1 = [t[-1] for t in taxa]  # Query consists of scientific names
    # t2 = [t[1] for t in taxa]  # Query consists of taxa ID
    # t3 = [t[0] for t in taxa]  # Query consists of OMA codes
    #
    # index = [1, 2] * 10
    # # Query consists of mix of scientific names and taxa IDs
    # t4 = [t[i] for t, i in zip(taxa, index)]
    #
    # index = [0, 1, 2] * 7
    # # Query consists of mix of scientific names, taxa IDs, and OMA codes
    # t5 = [t[i] for t, i in zip(taxa, index)]
    # oma(t5, outdir=folder)
