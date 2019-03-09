#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import logging
import argparse

from collections import namedtuple
from scipy.stats import poisson

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iTOL]')
warn, info, error = logger.warning, logger.info, logger.error


def _tester(obs, exp, values, alpha=0.05):
    """
    One sample T-test to determine whether the observed value is statistically
    significantly different to the expected value.

    :param obs: int, observed value.
    :param exp: float, the expected value.
    :param values: list or tuple, a list of expected values where exp was
        calculated.
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


def detect(branchpair=None, pars=None, cons=None, wd='', fn='', tester=None,
           printout=True, verbose=True):
    """
    Pairwise comparison for parallel and convergent amino acid replacements in
    protein sequences.

    :param branchpair: list, a list of branch pairs need to be tested.
    :param pars: dict, a dict object stores parallel changes.
    :param cons: dict, a dict object stores convergent changes.
    :param wd: str, path to the work directory. Without specifying, it will be
        set to current work directory. A file ends with '.counts.tsv' in the
        work directory will be used if neither pars nor cons was provided.
    :param fn: str: path to the result file.
    :param tester: function, a function for test the differences.
    :param printout: bool, print out the test results (default) or only return
        the test result without printing them out.
    :param verbose: bool, invoke verbose or silent process mode,
        default: False, silent mode.
    :return: list, a list of test results.
    """

    logger.setLevel(logging.INFO if verbose else logging.ERROR)
    
    if pars is None or cons is None:
        wd = wd if wd else os.getcwd()
        if wd and os.path.isdir(wd):
            fns = glob.glob(os.path.join(wd, '*.counts.tsv'))
            if fns:
                fn = fns[0]
            else:
                error('No result file was found, detect aborted.')
                sys.exit(1)
        elif fn:
            if not os.path.isfile(fn):
                error('Result {} is not a file or does not exist, detect '
                      'aborted.'.format(fn))
                sys.exit(1)
        else:
            error('The wd {} is not a directory or does not exist, detect '
                  'aborted.'.format(wd))
            sys.exit(1)

        pars, cons = {}, {}
        with open(fn) as f:
            for line in f:
                blocks = line.strip().split()
                if blocks[0] == 'P':
                    pars[blocks[1]] = blocks[2:]
                elif blocks[0] == 'C':
                    cons[blocks[1]] = blocks[2:]
                    
    elif not isinstance(pars, dict) and not isinstance(cons, dict):
        error('Invalid pars and cons, they need to be dict objects, detect '
              'aborted.')
        sys.exit(1)
    
    if isinstance(branchpair, str):
        pairs = [branchpair]
    elif isinstance(branchpair, (list, tuple)):
        pairs = branchpair
    else:
        info('No interested branch branchpair assigned, doing test for all '
             'branch pairs.')
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
                warn('Invalid branch branchpair {} was ignored'.format(p))
    
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
        h1 = ''.join(['|', ' ' * 34, '|', 'Parallelism'.center(21), '|',
                      'Convergence'.center(21), '|'])
        line1 = ''.join(['+', '-' * 34, '+', '------+', '------+', '-------+',
                         '------+', '------+', '-------+'])
        h2 = ''.join(['|', 'Branch Pair'.center(34), '|'] + [
            ' Obs. | Exp. |P-value|'] * 2)
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
    
    output = os.path.join(wd, 'pairwise.comparison.tsv')
    with open(output, 'w') as o:
        o.write('branchpair\tP_Obs\tP_Exp\tp1\tC_Obs\tC_Exp\tp2\n')
        o.writelines(
                '{}\t{}\t{:.4f}\t{:2.1E}\t{}\t{:.4f}\t{:2.1E}\n'.format(*result)
                for result in results)
    info('Successfully saved comparison results to {}'.format(output))
    return results


def main():
    des = """Pairwise comparison for parallel and convergent amino acid
replacements in protein sequences"""
    
    epilog = """
Only support Poisson test or one-sample t test, if you expect other tests,
please implement the test by yourself and using detect() in Python instead of
command line.
"""
    
    formatter = argparse.RawDescriptionHelpFormatter
    parse = argparse.ArgumentParser(description=des,
                                    prog='detect',
                                    usage='%(prog)s RESULT [OPTIONS]',
                                    formatter_class=formatter, epilog=epilog)
    parse.add_argument('RESULT',
                       help='Path to the result file contains identified '
                            'parallel and convergent changes.')
    parse.add_argument('-b', '--branchpair',
                       help='Comma separated branch pairs.')
    parse.add_argument('-v', '--verbose', action='store_true',
                       help='Invoke verbose or silent (default) process mode.')
    
    args = parse.parse_args()
    result, bp = args.RESULT, args.branchpair
    
    detect(bp, fn=result, verbose=args.verbose)
    

if __name__ == '__main__':
    main()
