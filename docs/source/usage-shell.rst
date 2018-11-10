.. _intro-usage-shell:


Using ProtParCon in terminal
============================

Check ProtParCon command-line toolsets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After ProtParCon has been installed, you should also have a six command-line 
tools installed, they are: msa, asr, mlt, aut, imc, and sim. To check the 
availability and usage of any of these command-line tools, you can type 
the name of the command-line tool with '-h' flag in a terminal to see the 
usage page, for example::

    $ msa -h 
    
The above command should print out the usage of `msa` without any error.
The '-h' flag can also be used to check the usage for other commands. We 
suggest you use this flag to display the usage of the command and learn 
how to use it from its usage.

Using toolsets in terminal
~~~~~~~~~~~~~~~~~~~~~~~~~~

By checking usage of command line tool, it should be very easy for you to use 
the toolsets shipped by ProtParCon. For example, the following command will 
align a sequence file named `seq.fa` using `MUSCLE` and save the alignment 
output to a file named `alignment.fasta`::

    $ msa muslce seq.fa -o alignment.fasta
    
In the following example, `imc` is used to automate sequence alignment, 
ancestral states reconstruction, sequence 
simulation, and identify parallel and convergent amino acid replacements 
both in the reconstructed ancestral sequences and the simulated sequences::

    $ imc seq.fa tree.newick -aligner muscle -ancestor codeml -simulator seqgen
    
After running the above command, there are several files stored in your current
work directory. You can get the information about identified parallel and 
convergent amino acid replacements in file `imc.counts.tsv` and file 
`imc.details.tsv`.

All command-line tools have nearly the same signatures as the equivalent
functions in python module, refer to the usage of each command and the 
examples showing the usage of the equivalent functions, it should be very
easy for you to use the command-line toolsets.

In the following part, we list examples of using command-line commands 
to do the same work that we have done in the ProtParCon python usage part. We are not
going to repeat the details of each example, if you have any question about
to details, see ProtParCon python usage part.

Align sequences using MUSCLE and save the alignment output to a file named 
`seq.muscle.fasta` (default name)::

    $ msa muscle seq.fa
    
Align sequence using MAFFT and save the alignment output to a FASTA format file 
named 'seq.mafft.fasta' (default name)::
    
    $ msa mafft seq.fa
    
And this will align the same sequence with Clustal (Omega) and save the
alignment to a FASTA file named 'seq.clustal.fasta'::

    $ msa clustal seq.fa

.. note::

    The above example assumes that you have system wide Clustal Omega
    installed and the string clustal point to the executable of Clustal
    Omega. If your Clustal Omega is not system widely installed, or the
    path to its executable is not `clustal`, change it accordingly.

Align sequence using MUSCLE and save the alignment output into a file name
`alignment.fasta`::

    msa muscle seq.fa -o alignment.fasta


Infer ML tree using IQ-TREE and save the best ML tree to a file named 
`msa.IQ-TREE.ML.newick` (default name)::

    $ mlt iqtree msa.fa
    
Infer ML tree using RAxML and save the ML tree to a file named
'msa.RAxML.ML.newick' (default name)::
    
    $ imc-mlt raxml msa.fa
    
Infer ML tree using PjyML and save the ML tree to a file named
'msa.PhyML.ML.newick' (default name)::
    
    $ mlt phyml msa.fa
    
Infer ML tree using FastTree and save the ML tree to a file named
'msa.FastTree.ML.newick' (default name)::
    
    $ mlt fasttree msa.fa
    
Infer ML tree using RAxML and save the tree into a file named `tree.newick`
in the same directory of the alignment file with '-o' option::

    $ mlt raxml msa.fa -o tree.newick


Infer ML tree using LG model with 8 Gamma categories accounting for among-site
rate variation and estimating ML base frequencies of 20 amino acids via PhyML::

    $ mlt PhyML msa.fa -o tree.newick -model LG+G8+F

The same as the above example, use '-gamma' and '-freq' options::

    $ mlt PhyML msa.fa -o tree.newick -model LG -gamma 9, -freq estimate

Infer ML tree with a start tree and/or constraint tree via `start_tree` and
`constraint_tree` options::

    $ mlt raxml msa.fa -o tree.newick -model LG+G8+I -stree /path/to/the/start/tree/file -ctree /path/to/the/constraint/tree/file


Reconstruct ancestral states using CODEML and save the ancestral states output
to a file named `msa.codeml.tsv` (default name)::

    $ asr codeml msa.fa tree.newick
    
Reconstruct ancestral states using RAxML and save the ancestral states
output to a file named 'msa.raxml.tsv' (default name)::
    
    $ asr raxml msa.fa

Reconstruct ancestral states using RAxML and save the ancestral states output
to a file named `ancestors.tsv` via '-o' option::

    $ asr raxml msa.fa -o ancestors.tsv
    
Reconstruct ancestral states via CODEML using WAG model (-model option)::

    $ asr codeml msa.fa tree.newick -model WAG

Reconstruct ancestral states via RAxML using 'WAG' model::

    $ asr raxml msa.fa tree.newick -model WAG
    
Reconstruct ancestral states using LG model with 8 Gamma categories and
a ML estimate of base frequencies of 20 amino acids via RAxML::

    $ asr raxml msa.fa tree.newick -model LG+G8+F
    
Do the same thing as the above example, but use '-gamma' and '-frq' options::
     
    $ asr raxml msa.fa tree.newick -model LG -gamma 8 -freq estimate

Use a specified model (or matrix) file along with complicated modeling
information for ancestral states reconstruction::
    
    $ asr codeml msa.fa tree.newick -model /path/to/my/own/model -gamma=8 -freq estimate

.. note::

    The model (or matrix) file needs to be in the right format required by ASR
    programs, before use the model file, check the manual for your ASR
    program to make sure you model file is in the right format.


Simulate sequences in the simplest way::

    $ sim evolver tree.newick

Use Seq-Gen to simulate 200 protein datasets with the length set to 500 amino
acids and substitution model set to LG with 8 Gamma categories to account for
among sites rate variation::

    $ sim seqgen tree.newick -length 500 -n 200 -model LG -gamma 8

Use Seq-Gen to simulate 200 protein datasets with the length and base
frequencies of 20 amino acids extracted from a multiple protein sequence
alignment file::

    $ sim seqgen tree.newick -n 200 -model LG -gamma=8 -msa /path/to/the/multiple/sequence/alignment/file -freq estimate

Topology test (AU test) using imc-aut::

    $ aut iqtree msa.fa tree.newick -model WAG


Identify parallel and convergent amino acid replacements using ancestral states
reconstruction generated by ProtParCon::

    $ imc path/to/the/ancestral/states/file
