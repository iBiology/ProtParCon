#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

VERSION = '1.0'

wd = os.path.abspath(os.path.dirname(__file__))


def readme():
    with open(os.path.join(wd, 'README.rst')) as f:
        return f.read()


setup(name='iMC',
      version=VERSION,
      description="""ProtParCon - A framework for framework for processing
      molecular data and identifying parallel and convergent amino acid
      replacements.""",
      long_description=readme(),
      url='https://github.com/iBiology/ProtParCon',
      author='FEI YUAN',
      author_email='yuanfeifuzzy@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['biopython>=1.71'],
      entry_points={
          'console_scripts': ['imc=imc.imc:main',
                              'oma=imc.oma:main',
                              'msa=imc.msa:main',
                              'mlt=imc.mlt:main',
                              'asr=imc.asr:main',
                              'aut=imc.aut:main',
                              'sim=imc.sim:main']
          },
      include_package_data=True,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          ],
      keywords='phylogeny tree alignment simulation biology bioinformatics'
      )


if __name__ == '__main__':
    pass
