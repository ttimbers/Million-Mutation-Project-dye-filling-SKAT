# Million-Mutation-Project-dye-filling-SKAT

Genomic data and code to accompany the SKAT analysis of Million Mutation Project strains 
dye-filling phenotypes reported in [Timbers *et al., bioRxiv*, 2015](http://dx.doi.org/10.1101/027540).
Specifically, this is the code used to generate Supplementary Tables 3-6.

All other data used for this project are stored in a private repository on Figshare, which will 
be made public upon publication of the manuscript.

## Dependencies

* Bash Shell (version 3.2.57(1))
* Make (version 3.81)
* Perl (v5.18.2)
* Plink (v1.90b1g)
* R (version 3.2.3)
* R packages: 
	1. dplyr (version 0.4.3)
	2. fdrtool (version 1.2.15)
	3. plyr (version 1.8.3) 
	4. pwr (version 1.1-3)
	5. SKAT (version 1.0.9) 
	6. stringr (version 1.0.0)  

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