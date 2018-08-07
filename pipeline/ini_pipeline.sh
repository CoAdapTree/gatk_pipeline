#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --job-name=initiatepipe
#SBATCH --export=all
#SBATCH --time=00:59:00
#SBATCH --mem=500mb
#SBATCH --cpus-per-task=1


python $HOME/pipeline/01a_trim-fastq.py $1 $2
