#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=8

# load docker image
udocker load docker_images/mmp-dyf-skat.tar

# make clean Million-Mutation-Project-dye-filling-SKAT in ttimbers/mmp-dyf-skat
udocker run --rm ttimbers/mmp-dyf-skat bash /home/ttimbers/projects/Million-Mutation-Project-dye-filling-SKAT/make_with_udocker.sh clean /home/ttimbers/projects/Million-Mutation-Project-dye-filling-SKAT