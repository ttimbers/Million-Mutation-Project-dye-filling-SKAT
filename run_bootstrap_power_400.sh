#!/bin/bash
#PBS -l walltime=00:20:00
#PBS -l nodes=1:ppn=8
#PBS -M ttimbers@sfu.ca

# load docker image
udocker load --input docker_images/mmp-dyf-skat.tar

# run Makefile to make all targets for Million-Mutation-Project-dye-filling-SKAT in ttimbers/mmp-dyf-skat
udocker run --rm -w /home/ttimbers/projects/Million-Mutation-Project-dye-filling-SKAT ttimbers/mmp-dyf-skat bash /home/ttimbers/projects/Million-Mutation-Project-dye-filling-SKAT/bin/bootstrap_power.sh '400'
