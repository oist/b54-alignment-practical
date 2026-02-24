#!/bin/bash
#SBATCH --job-name="iqtree"
#SBATCH --output="logs/%x_%A_%a.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-00:30:00
#SBATCH --array=0-9

module load bioinfo-ugrp-modules
module load Other/iqtree

ALIGNMENTS=("pax6_family.clustalw.fasta" "pax6_family.linsi.fasta"  "pax6_family.muscle.fasta"  "pax6_family.prank.fasta.best.fas" \
            "pax6_family.linsi.fasta.BMGE"  "pax6_family.muscle.fasta.BMGE"  "pax6_family.prank.fasta.best.fas.BMGE" \
            "pax6_family.linsi.fasta.trimal"  "pax6_family.muscle.fasta.trimal"  "pax6_family.prank.fasta.best.fas.trimal" )

ALIGNMENT_TO_WORK_WITH=${ALIGNMENTS[SLURM_ARRAY_TASK_ID]}

echo "Infering phylogeny of ${ALIGNMENT_TO_WORK_WITH}..."
time iqtree -T 4 -s $ALIGNMENT_TO_WORK_WITH -m LG+G

