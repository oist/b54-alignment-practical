#!/bin/bash
#SBATCH --job-name="diamond_drosophila_eyeless_isoformA_aa_seq_with_ncbinr"
#SBATCH --output="logs/%x_%j.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-2:00:00

module load bioinfo-ugrp-modules
module load Other/DIAMOND
module load DB/diamondDB/ncbi/2022-07

time diamond blastp --db ${DIAMONDDB}/nr.dmnd -q drosophila_pax6_isoformA.faa -p 16 -o result_drosophila_pax6_isoformA_diamond_ncbinr.out --header

