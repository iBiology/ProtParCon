Empirical model of protein evolution
====================================

This directory stores 12 empirical models of protein evolution adopted from `PAML 4.9e <http://abacus.gene.ucl.ac.uk/software/paml.html>`_:
    * cprev10
    * cprev64
    * dayhoff
    * dayhoff-dcmut
    * jones-dcmut
    * jtt
    * lg
    * mtart
    * mtmam
    * mtrev24
    * mtzoa
    * wag

All these models were save in plain text file without extension and the all have the same structure:
    * Line 1 is intentional leave as a blank line
    * Line 2 - 20 store symmetrical part (lower triangular) of the rate matrix (strictly enforcement)
    * Line 21 is intentional leave as a blank line (strictly enforcement)
    * Line 22 stores the amino acid frequencies (strictly enforcement)
    * Line 23 is intentional leave as a blank line (strictly enforcement)
    * Line 24 stores the order of the amino acid using one letter code (not necessary)
    * Line 25 stores the order of the amino acid using three letter code (not necessary)
    * Line 26 is a blank new line (not necessary)
    * For alignment purpose, data can be separated by any numbers of whitespace

All the unnecessary information about those models were deleted, see the original models in PAML 4.9e
package for more details.

All these models can be used in ``iparcon`` by their names (case insensitive) without specify their
relative or absolute pathname. If users try to use a model not listed here, however, the relative or
absolute pathname of the model file is needed.

If a model was used for ancestral states reconstruction, the original model file shipped with PAML 4.9e
or after should be fine, but a customised model file needs to be modified and make make sure the file
has the mentioned structure for symmetrical part (lower triangular) of the rate matrix and the amino
acid frequencies. But if the model file was used for simulation purpose, whether it is shipped by PAML
or not, users should always make necessary modification to make sure the file has required structure.

For Seq-Gen, it requires 190 rates representing the upper (off-diagonal) triangle of a 20x20 matrix
with amino acids in the order ARNDCQEGHILKMFPSTWYV, all models listed here cannot be directly used
in Seq-Gen, users need to convert the matrix first. After convert the lower triangle to upper triangle,
the 20 and 21 lines need to be blank lines and the amino acid frequencies need to be on the 22 line.
