.. _intro-usage-python:

The best way to learn is with examples, and ProtParCon is no exception. For this
reason, the following sections describe some pre-made examples of how to 
use ProtParCon to manipulate molecular sequence data and identify molecular 
convergences.

Since ProtParCon is written in pure Python, you can use ProtParCon in any Python script
or interactive Python session. Meanwhile, ProtParCon also ships a command-line
toolset to you, you are able to use the command-line tools in a terminal.
We will show you examples of using ProtParCon in Python and terminal separately.

Use ProtParCon in Python
========================

Import ProtParCon
~~~~~~~~~~~~~~~~~

You can import ProtParCon as a module or directly import general use functions
inside the module, thus the following two ways of importing are acceptable:

    >>> import ProtParCon

    >>> from ProtParCon import msa, aut, mlt, asr, imc

Align multiple sequences
~~~~~~~~~~~~~~~~~~~~~~~~

Multiple sequence alignment (MSA) can be easily done in ProtParCon by using
function ``msa()``, and this function has a common interface for all
supported MSA programs.

At this stage, MSA via MUSCLE, MAFFT, and Clustal (Omega) are supported.
You are only required to pass the path to the executable of any supported
MSA program and the path to a sequence file in FASTA format to function
``msa()``, ProtParCon will then take care of everything related to sequence
alignment for you.

The following examples assume that you have the needed MSA program installed
in your system and the path to its executable is the same as used here.
Otherwise, you need to replace the path to the executable of a MSA program
with its actual path. In our case, the executable of MUSCLE is named `muscle`
and it is already in the system path, so we are able to directly use `muscle`
as the path to the MUSCLE executable. All the examples also assume that the
sequence file used is in the current work directory and it is named `seq.fa`.
If this is not the case for you, you need to use either a relative or absolute
path of the files accordingly to make ProtParCon works as expected.
 
This Python session will automatically align sequences using MUSCLE and save
the alignment output to a file named `seq.muscle.fasta` in the same directory
of the sequence file:

    >>> from ProtParCon import msa
    >>> msa('muscle', 'seq.fa')
    
The default naming rule for alignment output will be in the format of
[basename].[aligner].fasta, where basename is the filename of the sequence
file without known FASTA format file's extension, aligner is the name of the
MSA program (figured by ProtParCon according to the MSA's executable) in lower
case, and fasta is the file extension for FASTA file. Thus, the following
session will align sequence using MAFFT and save the alignment output to a
FASTA format file named 'seq.mafft.fasta':
    
    >>> from ProtParCon import msa
    >>> msa('mafft', 'seq.fa')
    
And this will align the same sequence with Clustal (Omega) and save the
alignment to a FASTA file named 'seq.clustal.fasta':

    >>> from ProtParCon import msa
    >>> msa('clustal', 'seq.fa')

.. note::
    The above example assumes that you have system wide Clustal Omega
    installed and the string clustal point to the executable of Clustal
    Omega. If your Clustal Omega is not system widely installed, or the
    path to its executable is not `clustal`, change it accordingly.

If you want to save the alignment output to a file using a customized name,
you can pass a filename to argument ``outfile``. This session will align
sequence using MUSCLE and save the alignment output into a file named
`alignment.fasta` in the same directory of the sequence file (if you want
to save the alignment to a different directory, pass a relative or absolute
path accordingly):

    >>> from ProtParCon import msa
    >>> msa('muscle', 'seq.fa', outfile='alignment.fasta')


Infer Maximum-Likelihood tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Infer Maximum-Likelihood (ML) tree can be easily done in ProtParCon by using
function ``mlt()``. Since this function provides a common interface for
all supported ML programs, you can use this function to infer ML tree
using all supported ML programs in the same manner.

At this stage, ML tree inference via IQ-TREE, FastTree, RAxML, and PhyML
are supported in ProtParCon. You are only required to pass the path to the
executable of any supported MSA program and the path to an alignment file
in FASTA format to function ``msa()``, ProtParCon will then take care of
everything about ML tree inference for you.

The following examples assume you have the needed ML program installed on
your system and the path to its executable is the same as used here.
Otherwise, you need to replace the path to the executable of a MSA program
with its actual path. In our case, the executable of IQ-TREE is named
`iqtree` and it is already in the system path, so we are able to directly
use `iqtree` as the path to the executable of IQ-TREE software. All the
examples also assume that the alignment file used is in the current work
directory and it is named `msa.fa`. If this is not the case for you,
you need to use either a relative or absolute path of the files accordingly
to make ProtParCon works as expected.
 
This Python session will automatically infer ML tree using IQ-TREE and save
the best ML tree to a file named `msa.IQ-TREE.ML.newick` in the same
directory of the multiple sequence alignment file:

    >>> from ProtParCon import mlt
    >>> mlt('iqtree', 'msa.fa')
    
The default naming rule for ML tree output will be in the format of
[basename].[mlter].ML.newick, where basename is the filename of the
alignment file without known FASTA format file's extension, mlter is
the name of the ML program (figured by ProtParCon according to the ML program's
executable), ML denotes that the file contains a Maximum-Likelihood
tree, and newick is the file extension for NEWICK file. Thus, the following
session will infer ML tree using RAxML and save the ML tree to a file named
'msa.RAxML.ML.newick':
    
    >>> from ProtParCon import mlt
    >>> mlt('raxml', 'msa.fa')
    
Since ``mlt()`` provides a common interface for inferring ML tree, the
other two supported ML programs can also be used in the same way. This
example will infer ML tree using PjyML and save the ML tree to a file named
'msa.PhyML.ML.newick':
    
    >>> from ProtParCon import mlt
    >>> mlt('phyml', 'msa.fa')
    
And, of course, this example will infer ML tree using FastTree and save the
ML tree to a file named 'msa.FastTree.ML.newick':
    
    >>> from ProtParCon import mlt
    >>> mlt('fasttree', 'msa.fa')
    
If you want to save the resulted ML tree to a file using a customized
name, you can pass a filename to argument ``outfile``. This session will
infer ML tree using RAxML and save the tree to a file named `tree.newick`
in the same directory of the alignment file (if you want to save the
ancestral states output to a different directory, pass a relative or absolute
path accordingly):

    >>> from ProtParCon import mlt
    >>> mlt('raxml', 'msa.fa', outfile='tree.newick')


Function `mlt()` also allows you to specify a substitution model for inferring
a ML tree. The model can be a name of empirical substitution name (e.g. JTT,
WAG, LG, ...) or a name combined with a other options (e.g. JTT+G8,
LG+G8+I, WAG+G4+I+F). For example, the following example will infer ML tree
using LG model with 8 Gamma categories account for among-site rate variation
and estimated ML base frequencies of 20 amino acids via PhyML:

    >>> from ProtParCon import mlt
    >>> mlt('raxml', 'msa.fa', outfile='tree.newick', model='LG+G8+F')

If you do not want to specify the modeling process via argument `model`, there
are four other arguments: `gamma`, `alpha`, `freq`, and `invp`, that you can
use to specify additional modeling information, such as number of discrete
Gamma categories, shape of discrete Gamma distribution, base frequencies of
20 amino acids, and proportion of invariable site. The value specified by
these arguments have high priority than the ones coming with model argument.
Therefore, the following example will do the exactly same things as the above
example:

    >>> from ProtParCon import mlt
    >>> mlt('raxml', 'msa.fa', outfile='tree.newick', model='LG',
            gamma=9, freq='estimate')

You also have the option to provide a start tree and/or constraint tree to
control the way in which ML tree has been inferred via arguments `start_tree` and
`constraint_tree`:

    >>> from ProtParCon import mlt
    >>> mlt('raxml', 'msa.fa', outfile='tree.newick', model='LG+G8+I',
             start_tree='/path/to/the/start/tree/file',
             constraint_tree='/path/to/the/constraint/tree/file')



Reconstruct ancestral states
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ancestral states reconstruction (ASR) can be easily done in ProtParCon by using
function ``asr()``, and this function has a common interface for all supported
ASR programs.

At this stage, ASR via CODEML (inside PAML package), and RAxML are supported.
You are only required to pass the path to the executable of any supported ASR
program, the path to an alignment file in FASTA format, and a guide tree file
in NEWICK format to function ``asr()``, ProtParCon will then take care of everything
about ancestral states reconstruction and results parsing for you.

The following examples asuume you have the needed ASR program installed on your
system and the path to its executable is the same as used here. Otherwise, you 
need to replace the path to the executable of a ASR program with its actual
path. In our case, the executable of CODEML is named `codeml` and it is already
in the system path, so we are able to directly use `codeml` as the path to the
CODEML executable. All the examples also assume that the alignment file and the
tree file used here are in the current work directory and they are named
`msa.fa` and `tree.newick`. If this is not the case for you, you need to use
either a relative or absolute paths of these files accordingly to make ProtParCon
works as expected.
 
This Python session will automatically reconstruct ancestral states using CODEML
and save  the ancestral states output to a file named `msa.codeml.tsv` in the
same directory of the sequence file:

    >>> from ProtParCon import ars
    >>> asr('codeml', 'msa.fa', 'tree.nwick')
    
The resulted ancestral states file `seq.codeml.tsv` is TAB separated file, and
the first line of the file will start with '#TREE' and followed by a TAB (\t)
and then a NEWICK formatted tree string, the internal nodes are labeled. The
second line of the tsv file is intentionally left as a blank line and the rest
lines of the file are tab separated sequence IDs and amino acid sequences.
    
The default naming rule for ancestral states output will be in the format of
[basename].[asrer].tsv, where basename is the filename of the alignment file
without known FASTA format file's extension, asrer is the name of the ASR
program (figured by ProtParCon according to the MSA's executable) in lower case,
and fasta is the file extension for FASTA file. Thus, the following session
will reconstruct ancestral states using RAxML and save the ancestral states
output to a file named 'msa.raxml.tsv':
    
    >>> from ProtParCon import asr
    >>> asr('raxml', 'msa.fa')

If you want to save the ancestral states output to a file using a customized
name, you can pass the filename to argument ``outfile``. This session will
reconstruct ancestral states using RAxML and save the ancestral states output
to a file named `ancestors.tsv` in the same directory of the alignment file
(if you want to save the ancestral states output to a different directory,
pass a relative or absolute path accordingly):

    >>> from ProtParCon import asr
    >>> asr('raxml', 'msa.fa', outfile='ancestors.tsv')
    
The default reconstruction will use JTT model as substitution model, of course
you can use any supported substitution models by passing the name of the model
to argument ``model``. This session will reconstruct ancestral states via
CODEML using WAG model:

    >>> from ProtParCon import ars
    >>> asr('codeml', 'msa.fa', 'tree.nwick', model='WAG')

Since ``asr()`` function provides a common inteface for ancestral states
reconstruction, the same rule also applies to other support ASR programs,
like RAxML. Thus, this session will reconstruct ancestral states via RAxML
using 'WAG' model:

    >>> from ProtParCon import ars
    >>> asr('codeml', 'msa.fa', 'tree.nwick', model='WAG')
    
The argument `model` can also pass other information for guiding state
reconstruction rather than just the name of a empirical substitution model.
If the model argument combined name and Gamma category numbers, i.e. JTT+G4,
WAG+G8, etc., a discrete Gamma model would be used to account for among-site
rate variation. If the model argument combined name and frequencies, i.e.
LG+G8+F, WAG+F, etc., a ML estimate of base freqeuencies of 20 amino acids
would be used instead of empirical values associated with the specified
substitution model. Thus, this session will reconstruct ancestral states
using LG model with 8 Gamma categories and a ML estimate of base frequencies
of 20 amino acids via RAxML:

    >>> from ProtParCon import ars
    >>> asr('raxml', 'msa.fa', 'tree.nwick', model='LG+G8+F')
    
Another way for you to providing complicted modeling information is to use
these four arguments: `gamma`, `alpha`, `freq`, and `invp`, along with
argument `model`. Argument `gamma` will accept a integer of number of Gamma
category, `alpha` will let you specify the shape parameter of the discrete
Gamma distribution, and `freq` will accept either `empirical`
or `estimate` to specify how the base frequencies of 20 amino acids will be
handled. Thus, the following example will do the exactly same thing as the
above example:
     
    >>> from ProtParCon import ars
    >>> asr('raxml', 'msa.fa', 'tree.nwick', model='LG',
            gamma=8, freq='estimate')

Moreover, ProtParCon also allows you to use a customized substitution model (or
matrix) instead of the built-in empirical models in ASR programs. You can
pass the model (or matrix) file to argument `model`. In this case, if you
still want to provide modeling information, such as Gamma categories and
shape, base frequencies of amino acid, you are required to pass all these
information through `gamma`, `alpha`, and `freq`. This example shows you
how you can use a specified model (or matrix) file along with complicated
modeling information for ancestral states reconstruction:
    
    >>> from ProtParCon import ars
    >>> asr('codeml', 'msa.fa', 'tree.nwick', model='/path/to/my/own/model', 
            gamma=8, freq='estimate')

.. note::
    The model (or matrix) file needs to be in the right format required by ASR
    programs, before use the model file, check the manual for your ASR
    program to make sure you model file is in the right format.


Simulate protein sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~
ProtParCon supports simulate proteins sequences via `EVOLVER` (inside PAML package)
and Seg-Gen in Python script or interactive session using a common interface.
You are only required to provide an executable of any supported simulation
program and a phylogenetic tree (with branch lengths) to the general use
function ``sim()``, ProtParCon will then take care everything about simulation for
you. The simplest example will look like this:

    >>> from ProtParCon import sim
    >>> sim('evolver', 'tree.nwick')

The above code will use EVOLVER to simulate proteins sequences along the
phylogenetic tree specified in tree file `tree.newick`. ProtParCon will
automatically extract the number of leaf nodes that need to be simulated from the
tree file. The default simulation procedure will set the substitution model
to JTT model, the length of protein sequence to 100 amino acid sites, and the
number of datasets (or duplicates) to 100. Thus, after run the above code,
ProtParCon will simulate protein datasets and save simulated sequences to a tab
separated file. The default naming rule for simulation output will be in the
format of [simulator].simulations.tsv, where simulator is the name of the
simulation program you used (this will be figured by ProtParCon automatically).

The resulted simulation file is a TAB separated text file, and the first line
of the file will start with '#TREE' and followed by a TAB (\t) and then a
NEWICK formatted tree string, the internal nodes are labeled. The second line
of the tsv file is intentionally left as a blank line and the rest
lines of the file are blank line separated sequence blocks, in each sequence
block, each sequence takes a whole line and the sequence IDs and amino acid
sequences are separated by a TAB ('\t').


The function `sim()` also allows you to simulate sequences under various
scenarios. For example, the following example will use Seq-Gen to simulate
200 protein datasets with the length set to 500 amino acids and substitution
model set to LG with 8 Gamma categories to account for among sites rate
variation:

    >>> from ProtParCon import sim
    >>> sim('seqgen', 'tree.newick', length=500, n=200,
            model='LG', gamma=8)

The following example will use Seq-Gen to simulate
200 protein datasets with the length and base frequencies of 20 amino acids
extracted from a multiple protein sequence alignment file:

    >>> from ProtParCon import sim
    >>> sim('seqgen', 'tree.newick', n=200, model='LG', gamma=8,
            msa='path/to/the/multiple/sequence/alignment/file',
            freq='estimate',
            outfile='path/to/the/simulation/output/file')

Since you also passed a filename to argument `outfile`, in the above example,
simulated sequences and the labeled tree will be saved to the file you
specified.

Topology test
~~~~~~~~~~~~~

For our own purpose, the only topology test supported in ProtParCon at this stage
is AU test and the test need to use IQ-TREE. The test can be easily done via
function `aut` inside ProtParCon. For example, assume we have multiple sequence
alignment file and a hypothesis phylogenetic relationship specified in a
NEWICK tree file, we want to know how big the differences are between the
phylogenetic relationship specified by the tree file and the relationship
specified by the best ML tree given the MSA file, thus we can do the following:

    >>> from ProtParCon import aut
    >>> aut('iqtree', 'msa.fa', 'tree.newick', model='WAG')

Running the above code will conduct AU test implemented in IQ-TREE program
without knowing how to use IQ-TREE.


Identify molecular convergence in proteins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function ``imc()`` inside ProtParCon package provides a common and easy way to
identify parallel and convergent amino acid replacements in orthologous
proteins separately. It is able to take the following kinds of
sequence data and identify parallel and convergent amino acid replacements
for all comparable branch pairs within the phylogenetic tree:

    * sequences: raw protein sequence file, need to be in FASTA format
      (a NEWICK format tree, a multiple sequence alignment program,
      an ancestral states reconstruction program are also required.
    * msa: multiple sequence alignment file, need to be in FASTA format
      (a NEWICK format tree and an ancestral states reconstruction program
      are also required).
    * ancestors: reconstructed ancestral states file, need to be in tsv
      (tab separated) file, the first line needs to start with #TREE,
      second line need to be a blank line, and the rest lines in the
      file need to be tab separated sequence name (or ID) and amino
      acid sequences.
    * simulations: simulated sequences, need to be in tsv file, the
      first line needs to start with #TREE, second line need to be
      a blank line, each dataset need to be separated by a blank line
      and inside each dataset block, each line should consist of tab
      separated sequence name (or ID) and amino acid sequences.


In the simplest way, if you pass a result file generated during ancestral
states reconstruction to `imc()`, parallel and convergent amino acid
replacements will be identified without requiring any other information:

    >>> from ProtParCon import imc
    >>> imc('path/to/the/ancestral/states/file')

After running the above code, there are two TAB separated files saved
to the current work directory: 'imc.counts.tsv' and 'imc.details.tsv'. The
former file contains information about the number of parallel and convergent
amino acid replacements among branch pairs, while the later file contains
details of each identified parallel or convergent amino acid replacement.

If you are interested in identifying parallel and convergent amino acid
replacements in a protein orthologous group and you also want to identify
the corresponding number of parallel changes in simulated datasets
based on the protein orthologous group, you can manually do sequence
alignment, ancestral states reconstruction, sequence simulation, and identify
parallel and convergent changes step by step, or you can just let `imc()`
do all the work for you:

    >>> from ProtParCon import imc
    >>> imc('path/to/the/orthologous/file', tree='path/to/the/tree/file',
            aligner='path/to/the/executable/of/a/MSA/program',
            ancestor='path/to/the/executable/of/a/ASR/program',
            anc_model='LG+G8+F',
            simulator='path/to/the/executable/of/a/simulation/program',
            sim_model='LG+G8+F',
            wd='path/to/the/work/directory')


