#!/bin/bash
#SBATCH -J  PauliShuffle     # Job name    
#SBATCH -N  1                  # Total Nodes 
#SBATCH -n  1                  # Total Tasks = Total Nodes
#SBATCH -p  gpu              # Queue: maxnodes = 32
#SBATCH -o  output/PauliShuffle%j.o # Stdout file
#SBATCH -e  output/PauliShuffle%j.e # Stderr file
#SBATCH -t  02:00:00         # Run time (hh:mm:ss)  

# launch the launcher
module load launcher
export TACC_LAUNCHER_PPN=1
export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
export WORKDIR=./output/$SLURM_JOB_ID/
export CONTROL_FILE=paramlist
export TACC_LAUNCHER_SCHED=interleaved 

$TACC_LAUNCHER_DIR/paramrun SLURM $EXECUTABLE $WORKDIR $CONTROL_FILE

echo "Job Complete"
