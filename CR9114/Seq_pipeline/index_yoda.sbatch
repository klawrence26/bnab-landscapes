#! /bin/bash -e
#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J yoda_index # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 4000 # Memory in MB
#SBATCH -o slurm_outputs/index_%A_%a.out
#SBATCH -e index_%A_%a.err

CURRDIR=$PWD


grep -P "1\t1\t1\t1\t" parsed_files/geno_9114_a_7_1 | sort  > corrected_files/index_9114_a_7_1.txt
grep -P "1\t2\t1\t2\t" parsed_files/geno_9114_a_7_2 | sort  > corrected_files/index_9114_a_7_2.txt
grep -P "1\t3\t1\t5\t" parsed_files/geno_9114_a_7_3 | sort  > corrected_files/index_9114_a_7_3.txt
grep -P "1\t4\t1\t6\t" parsed_files/geno_9114_a_7_4 | sort  > corrected_files/index_9114_a_7_4.txt

