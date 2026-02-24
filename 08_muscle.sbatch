#!/bin/bash
#SBATCH --job-name="muscle"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-2:00:00

module load bioinfo-ugrp-modules
module load DebianMed
module load muscle

time muscle -align pax6_family.fasta -output pax6_family.muscle.fasta

