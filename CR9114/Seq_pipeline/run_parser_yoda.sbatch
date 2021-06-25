#! /bin/bash -e
#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J parse # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-20:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 4000 # Memory in MB
#SBATCH -o slurm_outputs/parse_%A_%a.out
#SBATCH -e slurm_outputs/parse_%A_%a.err

CURRDIR=$PWD

bash clean_YODA fastq/9114/9114_$1_$2_${SLURM_ARRAY_TASK_ID}_S*_L001_R1_001.fastq.gz fastq/9114/9114_$1_$2_${SLURM_ARRAY_TASK_ID}_S*_L001_R2_001.fastq.gz | python yoda_parser.py > parsed_files/9114_$1_$2_${SLURM_ARRAY_TASK_ID}_L1
bash clean_YODA fastq/9114/9114_$1_$2_${SLURM_ARRAY_TASK_ID}_S*_L002_R1_001.fastq.gz fastq/9114/9114_$1_$2_${SLURM_ARRAY_TASK_ID}_S*_L002_R2_001.fastq.gz | python yoda_parser.py > parsed_files/9114_$1_$2_${SLURM_ARRAY_TASK_ID}_L2


