#!/bin/bash

rm *.out
rm *.b_dat
NUM_OMP=1 NUM_MPI=1 ./createJobscript.sh
sbatch job.script
