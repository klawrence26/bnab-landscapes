#! /bin/bash -e
#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J yoda_nuc_to_aa # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-02:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 4000 # Memory in MB
#SBATCH -o slurm_outputs/nuc_to_aa_%A_%a.out
#SBATCH -e slurm_outputs/nuc_to_aa_%A_%a.err

CURRDIR=$PWD

cat parsed_files/9114_$1_$2_${SLURM_ARRAY_TASK_ID} | python nuc_to_aa.py bin | sed '/N/d' > parsed_files/geno_9114_$1_$2_${SLURM_ARRAY_TASK_ID}

