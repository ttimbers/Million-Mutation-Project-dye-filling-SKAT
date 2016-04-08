# !bin/bash
# Tiffany Timbers, April 8, 2016
#
# Bash script to run with udocker (or docker) & ttimbers/mmp-dyf-skat. It functions to add 
# plink to $PATH and run Makefile in Million-Mutation-Project-dye-filling-SKAT directory
#
# Arguments: 1. path to Million-Mutation-Project-dye-filling-SKAT directory inside container
#
# How to run on Compute Canada's guillimin server:
# udocker run --rm ttimbers/mmp-dyf-skat bash /home/Million-Mutation-Project-dye-filling-SKAT/make_with_udocker.sh /home/Million-Mutation-Project-dye-filling-SKAT
#
# Example of how I ran it on my local machine:
# docker run --rm -v /Users/tiffanytimbers/Documents/Post-Doc/Manuscripts/MMP_dyf_screen/code/Million-Mutation-Project-dye-filling-SKAT:/home/Million-Mutation-Project-dye-filling-SKAT ttimbers/mmp-dyf-skat bash /home/Million-Mutation-Project-dye-filling-SKAT/make_with_udocker.sh /home/Million-Mutation-Project-dye-filling-SKAT

# export path to plink to $PATH
export PATH="/usr/bin/plink:$PATH"

# run Makefile in Million-Mutation-Project-dye-filling-SKAT directory
make -C $1