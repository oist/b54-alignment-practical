#!/bin/bash
#SBATCH --job-name="mafft_linsi"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-2:00:00

module load bioinfo-ugrp-modules
# do not load DebianMed, because its mafft version is broken
module load mafft

time linsi --thread 16 pax6_family.fasta > pax6_family.linsi.fasta

