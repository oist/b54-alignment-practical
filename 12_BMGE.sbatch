#!/bin/bash
#SBATCH --job-name="BMGE"
#SBATCH --output="logs/%x_%A_%a.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=compute
#SBATCH --time=0-00:30:00
#SBATCH --array=0-3

module load java-jdk/11

ALIGNMENTS=("pax6_family.clustalw.fasta" "pax6_family.linsi.fasta"  "pax6_family.muscle.fasta"  "pax6_family.prank.fasta.best.fas")

ALIGNMENT_TO_WORK_WITH=${ALIGNMENTS[SLURM_ARRAY_TASK_ID]}

echo "Trimming ${ALIGNMENT_TO_WORK_WITH}..."
time java -jar BMGE/src/BMGE.jar -i $ALIGNMENT_TO_WORK_WITH -t AA -o ${ALIGNMENT_TO_WORK_WITH}.BMGE -oh ${ALIGNMENT_TO_WORK_WITH}.BMGE.html

