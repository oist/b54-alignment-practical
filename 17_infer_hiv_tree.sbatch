#!/bin/bash
#SBATCH --job-name="iqtree_hiv"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-01:30:00

module load bioinfo-ugrp-modules
module load Other/iqtree

time iqtree -T 4 -s hiv_samples.mafft.fasta -m MFP -B 1000

