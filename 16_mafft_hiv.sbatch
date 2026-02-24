#!/bin/bash
#SBATCH --job-name="mafft_auto_hiv"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-00:10:00

module load bioinfo-ugrp-modules
# do not load DebianMed, because its mafft version is broken
module load mafft

time mafft --thread 4 --auto hiv_samples.fasta > hiv_samples.mafft.fasta

