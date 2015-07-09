#!/bin/bash

rm -f *.out
rm -f *.b_dat
NUM_OMP=8 NUM_MPI=4 ./createJobscript.sh
sbatch job.script
