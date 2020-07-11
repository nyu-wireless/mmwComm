#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ps3857@nyu.edu

module purge
module load matlab/2020a

RUNDIR=$SCRATCH/mmwComm/projects/dl_throughput
cd $RUNDIR

export SLURM_ARRAY_TASK_ID=2

matlab -nodisplay -nodesktop -r "run rffeTest.m"
