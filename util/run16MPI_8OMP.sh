#!/bin/bash

rm *.out
rm *.b_dat
NUM_OMP=8 NUM_MPI=16 ./createJobscript.sh
sbatch job.script
