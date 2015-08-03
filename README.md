# Million-Mutation-Project-dye-filling-SKAT

Genomic data and code to accompany the SKAT analysis of Million Mutation Project strains 
dye-filling phenotypes reported in Timbers et al.

All data used for this project is stored in a private repository on Figshare, which will 
be made public upon publication of the manuscript.

## Dependencies

`Bash Shell`, `Make`, `Perl`, `Plink v1.90b1g`, `R` and `R packages SKAT (version 0.95), stringr, fdrtool and plyr`

This code assumes that you have put `Plink` in your `Bash Shell`'s `$PATH`. This can be 
done by adding the following to `.bash_profile` in your home directory (`cd ~/` to get 
there).

~~~
# added Plink
export PATH="/path/to/where/you/installed/plink:$PATH"
~~~

## How to use it

1. Clone this repository
2. Navigate the Shell to the root directory of this project
3. Type `make`

#### This code currently works, but the phenotype data has not yet been released, and so you will hit an error without providing this data. This will be done very shortly. 