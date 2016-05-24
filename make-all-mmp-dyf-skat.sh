#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=8

# load docker image
udocker load --input /home/ttimbers/docker_images/mmp-dyf-skat.tar

# run Makefile to make all targets for Million-Mutation-Project-dye-filling-SKAT in ttimbers/mmp-dyf-skat
#udocker run --rm -w /gs/project/qkh-103-aa/Million-Mutation-Project-dye-filling-SKAT ttimbers/mmp-dyf-skat bash /gs/project/qkh-103-aa/Million-Mutation-Project-dye-filling-SKAT/make_with_udocker.sh all
udocker run --rm -it -v /gs/project/qkh-103-aa/Million-Mutation-Project-dye-filling-SKAT:/home/Million-Mutation-Project-dye-filling-SKAT -w /home/Million-Mutation-Project-dye-filling-SKAT ttimbers/mmp-dyf bash /gs/project/qkh-103-aa/Million-Mutation-Project-dye-filling-SKAT/make_with_udocker.sh all
