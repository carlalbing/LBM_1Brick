#!/bin/bash

rm -f *.out
rm -f *.b_dat
NUM_OMP=1 NUM_MPI=1 ./createJobscript.sh
sbatch job.script
