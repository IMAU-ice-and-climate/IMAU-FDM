#!/bin/bash

#SBATCH --job-name=maxFib_2D      ## Name of the job
#SBATCH --output=maxFib_2D.out    ## Output file
#SBATCH --time=48:00:00           ## Job Duration
#BATCH --ntasks=1            ## Number of tasks (analyses) to run
#SBATCH --cpus-per-task=2      ## The number of threads the code will use
#SBATCH --mem-per-cpu=19999M     ## Real memory(MB) per CPU required by the job.

## Load the python interpreter
module load python3

## Execute the python script
srun python3 FGR055_make_2Ddetail_ERA5_FDMv12A.py
