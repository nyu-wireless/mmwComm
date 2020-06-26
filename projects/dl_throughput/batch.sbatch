#!/bin/bash

module purge
module load matlab/2020a

RUNDIR=$SCRATCH/mmwComm/projects/dl_throughput
cd $RUNDIR

matlab -nodisplay -nodesktop -r "run rffeTest.m"