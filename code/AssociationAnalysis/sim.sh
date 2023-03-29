#!/bin/bash
#SBATCH --job-name=file_name
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user-email
#SBATCH -t 48:0:00
#SBATCH --output=myjobarray_%A-%a.out
#SBATCH --array=1-500
#SBATCH --mem=16000
#SBATCH -p node
echo "$SLURM_ARRAY_TASK_ID"

module load R
Rscript file_name

