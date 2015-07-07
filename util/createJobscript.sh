#!/bin/bash

jobFile="job.script"

cpusPerNode=8

tasksPerNode=$((cpusPerNode/NUM_OMP))
check=$((tasksPerNode*NUM_OMP))
if [ $check != $cpusPerNode ]; then
  echo "Wrong NUM_OMP: $NUM_OMP -> $check"
  exit 1
fi

if [ $NUM_MPI -lt $tasksPerNode ]; then
  nodes=1
  tasksPerNode=$NUM_MPI
else
  nodes=$((NUM_MPI/tasksPerNode))
  check=$((nodes*tasksPerNode))
  if [ $check != $NUM_MPI ]; then
    echo "Wrong NUM_MPI: $NUM_MPI -> $check"
    exit 1
  fi
fi

echo "Creating jobscript for ${NUM_MPI} MPI procs with ${NUM_OMP} OMP threads:"
echo "    Tasks per node: $tasksPerNode"
echo "             Nodes: $nodes"


{   # create the job script
    echo "#!/bin/bash -l" 
    echo "#SBATCH --nodes=$nodes" 
    echo "#SBATCH --ntasks-per-node=$tasksPerNode" 
    echo "#SBATCH --cpus-per-task=$NUM_OMP" 
    echo "#SBATCH --time=00:30:00" 
    echo "#SBATCH --res=eurohack15" 
    echo "export OMP_NUM_THREADS=$NUM_OMP" 
    echo "aprun -B ./WMBrick3D" 
} > $jobFile
