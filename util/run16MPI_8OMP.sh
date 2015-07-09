#!/bin/bash

rm -f *.out
rm -f *.b_dat
NUM_OMP=8 NUM_MPI=16 ./createJobscript.sh
sbatch job.script
