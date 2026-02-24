#!/bin/bash
#SBATCH --job-name="blastn_drosophila_eyeless_na_seq_with_nt"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-1:00:00

module load ncbi-blast
module load bioinfo-ugrp-modules
module load DB/blastDB/ncbi/2025-12-07-nt

time blastn -db nt -query drosophila_pax6.fna -num_threads 8 -out result_drosophila_pax6_blastn.out

