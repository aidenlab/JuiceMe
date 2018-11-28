Scripts for reproducing Hi-Culfite plots.

# General usage #
Follow the script make_oe.sh to replicate the Hi-Culfite analysis shown in 
the paper. To create the co-methylation plot shown in the supplement, look 
at make_cometh_oe.py

# Neighborhood Observed and Expected #
For some a given resolution, or binning, let *M* be the matrix of contacts 
in which both read ends are methylated; let *U* be the matrix of contacts in
which both read ends are unmethylated; and let *Y* be the matrix of contacts in
which one read end is methylated and one read end is unmethylated. (Note that
*Y* is not symmetric.)

Let *a* be the methylation vector, defined below.

We want to create a matrix that tells us, for every locus i, the likelihood
that i is methylated conditioned on the loci it is in contact with. We call
this the positional observed. The expected matrix is simply a column vector
corresponding to the one-dimensional methylation vector, repeated for every column.

*O<sub>i,j</sub>=M<sub>i,j</sub>+Y<sub>i,j</sub>*

*E<sub>i,j</sub>=a(i)* for all *j*

*T<sub>i,j</sub>=M<sub>i,j</sub>+U<sub>i,j</sub>+Y<sub>i,j</sub>+Y<sub>j,i<sub>*

*a* is the methylation vector. It is defined as the sum of methylated reads in 
the bin divided by the total number of reads interacting with that bin.

*a(i) = (&Sigma;<sub>j</sub> (M<sub>i,j</sub> + Y<sub>i,j</sub>))/(&Sigma;<sub>j</sub> (M<sub>i,j</sub> + Y<sub>i,j</sub> + Y<sub>j,i</sub> + U<sub>i,j</sub>))*


# Methylation Correlation #
For the methylation correlation analysis, we would like to determine if the methylation state of a read correlates with the methylation state of the neighboring sequence. We define the methylation correlation as the frequency with which locus j is methylated given that locus i is methylated, divided by the total number of times locus j is methylated. This is 

*M<sub>i,j</sub> / (M<sub>i,j</sub>+Y<sub>j,i</sub>)*

Similarly, we examine the unmethylation correlation: the number of times j is unmethylated given that i is unmethylated, divided by the total number of times locus j is unmethylated. This is

*U<sub>i,j</sub>/ (U<sub>i,j</sub>+Y<sub>i,j</sub>) *


# Comethylation Observed and Expected #
Using the same definitions as above, we want to determine if contacts are 
co-methylated at a greater frequency than one would expect from a null model.
In this case, the observed is

*O<sub>i,j</sub>=(M<sub>i,j</sub>+U<sub>i,j</sub>) / T<sub>i,j</sub>*

The expected is

*E<sub>i,j</sub>=((a<sub>i</sub>&ast;a<sub>j</sub>)+((1-a<sub>i</sub>)&ast;(1-a<sub>j</sub>)))&ast;T<sub>i,j</sub>*

# Correlation Plots #
*methylKitAnalysis.txt* - R script for creating the correlation plots 
between WGBS and Hi-Culfite methylation tracks. The input files to this 
function are created by aligning reads via bwameth, marking duplicates that
were called by Juicer via mark_dups.awk, and running 
`MethylDackel -F1024 extract --methylKit <reference sequence> <bam>`


