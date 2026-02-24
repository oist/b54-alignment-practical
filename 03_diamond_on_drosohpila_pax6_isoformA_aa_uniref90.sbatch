#!/bin/bash
#SBATCH --job-name="diamond_drosophila_eyeless_isoformA_aa_seq_with_uniref90"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-2:00:00

module load bioinfo-ugrp-modules
module load Other/DIAMOND
module load DB/diamondDB/uniref90/2020_06

time diamond blastp --db ${DIAMONDDB}/uniref90.dmnd -q drosophila_pax6_isoformA.faa -p 16 -o result_drosophila_pax6_isoformA_diamond_uniref90.out --header

