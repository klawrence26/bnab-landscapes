#!/bin/bash
#
#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J yoda_counts # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 2000 # Memory in MB
#SBATCH -o slurm_outputs/counts_%A_%a.out
#SBATCH -e slurm_outputs/counts_%A_%a.err

CURRDIR=$PWD


cut -f 1 corrected_files/final_9114_$1_$2_${SLURM_ARRAY_TASK_ID}.txt | sort | uniq -c | sort -nr | awk "{print \$2,\$3,\$1}" OFS="\\t"  > count_files/counts_9114_$1_$2_${SLURM_ARRAY_TASK_ID}.txt
