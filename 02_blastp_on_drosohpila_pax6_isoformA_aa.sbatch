#!/bin/bash
#SBATCH --job-name="blastp_drosophila_eyeless_isoformA_aa_seq_with_uniref"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-1:00:00

module load ncbi-blast
module load bioinfo-ugrp-modules
module load DB/blastDB/uniref90/2020_06

time blastp -db uniref90 -query drosophila_pax6_isoformA.faa -num_threads 8 -out result_drosophila_pax6_isoformA_blastp_uniref90.out

