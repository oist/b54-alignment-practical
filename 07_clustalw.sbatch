#!/bin/bash
#SBATCH --job-name="clustalw"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-2:00:00

module load bioinfo-ugrp-modules
module load DebianMed
module load clustalw

time clustalw -infile=pax6_family.fasta -outfile=pax6_family.clustalw.fasta

