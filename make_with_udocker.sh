# !bin/bash
# Tiffany Timbers, April 8, 2016
#
# Bash script to run with udocker (or docker) & ttimbers/mmp-dyf-skat. It functions to add 
# plink to $PATH and run Makefile in Million-Mutation-Project-dye-filling-SKAT directory
#
# Arguments: 	1. target to make 
#				2. path to Million-Mutation-Project-dye-filling-SKAT directory inside container
#
# How to run on Compute Canada's guillimin server to make all:
# udocker run --rm ttimbers/mmp-dyf-skat bash /home/ttimbers/projects/Million-Mutation-Project-dye-filling-SKAT/make_with_udocker.sh all /home/ttimbers/projects/Million-Mutation-Project-dye-filling-SKAT
#
# Example of how I ran it on my local machine to make all:
# docker run --rm -v /Users/tiffanytimbers/Documents/Post-Doc/Manuscripts/MMP_dyf_screen/code/Million-Mutation-Project-dye-filling-SKAT:/home/Million-Mutation-Project-dye-filling-SKAT ttimbers/mmp-dyf-skat bash /home/Million-Mutation-Project-dye-filling-SKAT/make_with_udocker.sh all /home/Million-Mutation-Project-dye-filling-SKAT

# export path to plink to $PATH
export PATH="/usr/bin/plink:$PATH"

# run Makefile in Million-Mutation-Project-dye-filling-SKAT directory
make $1 -C $2