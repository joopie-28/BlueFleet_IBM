#!/bin/bash
#SBATCH --job-name=ScreenParams
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=qbiol
#SBATCH --mem=2GB              # memory (MB)
#SBATCH --time=1-00:00          # time (D-HH:MM)
#SBATCH -o MyJob.%N.%j.out     # STDOUT
#SBATCH -e MyJob.%N.%j.err     # STDERR

#SBATCH --array=1-3

echo "Start time: "; date

module load applications/R/4.0.3

Rscript --vanilla Hunting_Behaviour_Run.R $SLURM_ARRAY_TASK_ID

echo "End time: "; date
