#!/bin/bash
#SBATCH --job-name=ASTRAL
#SBATCH --mem=250G
#SBATCH -n 1 #tasks
#SBATCH -N 1 #nodes
#SBATCH -c 16 #cores per task
#SBATCH -o slurm_output-ASTRAL.%J
#SBATCH -e slurm_error-ASTRAL.%J
#SBATCH -p threaded
#SBATCH --qos threaded
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kmkocot@ua.edu

module load compilers/gcc/6.5.0 
module load java/1.8.0
module load gccutil

java -XX:ConcGCThreads=3 -Xmx250000M -jar -D"java.library.path=/kmk/scripts/ASTRAL/Astral/lib" /kmk/scripts/ASTRAL/Astral/astral.5.7.5.jar -b bootstrap_file_list.txt -i trees.tre -o output_species_tree.tre
