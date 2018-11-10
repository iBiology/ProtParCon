.. _intro-install:


Install ProtParCon
==================

ProtParCon runs only on Python 3.4 or above, if you are already
familiar with installation of Python packages, you can easily install
ProtParCon and its dependencies from PyPI with the following command on
all major platforms (Windows, MacOS, or Linux):

    pip install ProtParCon

It is strongly recommended that you install ProtParCon in a dedicated 
virtualenv, to avoid conflicting with your system packages.


Install ProtParCon to a virtual environment (recommended)
=========================================================

We recommend installing ProtParCon inside a virtual environment on all 
platforms.

Python packages can be installed either globally (a.k.a system wide), or in
a user specified space. We do not recommend installing ProtParCon system wide.

Instead, we recommend that you install ProtParCon within a so-called "virtual
environment" (`virtualenv`_). Virtualenvs allows users not to conflict with
already-installed Python system packages (which could break some of your
system tools and scripts), and still install packages normally with ``pip``
(without ``sudo`` and the likes).

To get started with virtual environments, see
`virtualenv installation instructions`_.

.. _virtualenv: https://virtualenv.pypa.io
.. _virtualenv installation instructions: https://virtualenv.pypa.io/en/stable/installation/


Things that are good to know
============================

ProtParCon is written in pure Python and depends only on standard libraries
and one third party library: `Biopython`_.

.. note::

    `NumPy`_ is required by biopython, however, if biopython can be
    successfully installed as a dependency package of ProtParCon, NumPy
    should not be a problem.

Although ProtParCon itself is platform independent, the built-in supports for
many programs, e.g. MAFFT for sequence alignment, FastTree for phylogenetic
tree inference, and Seq-Gen for sequence simulation may not be platform
independent, it is strongly recommended that ProtParCon is used on MacOS or 
Linux platforms.

ProtParCon is designed for providing a easy and common interface for various
programs related to phylogenetic analysis and analysis of molecular parallel
and convergent amino acid replacements, therefore, users are recommended to
use ProtParCon inside their scripts for efficiently calling external programs
and building their pipeline by chaining several external programs.

.. _Biopython: https://biopython.org
.. _NumPy: www.numpy.org/
