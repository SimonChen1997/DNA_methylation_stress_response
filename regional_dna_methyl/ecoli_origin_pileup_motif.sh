#!/bin/bash -l
#SBATCH --job-name="ecoli_origin_pileup_motif"
#SBATCH --partition=general
#SBATCH --array=25-60
#SBATCH -o ecoli_origin_pileup_motif.o
#SBATCH -e ecoli_origin_pileup_motif.e

#########################################################
module load anaconda3
module load samtools

#########################################################
### get modification at base level
source activate ont-modkit-0.5.0

########## e.coli
## targeting m6a
modkit pileup $minimap2_primary_m6a/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam $ecoli_motif_focus_pileup/barcode${SLURM_ARRAY_TASK_ID}_m6a_gatc.bed --motif GATC 1 --ref $ecoli_ref

## targeting m4c
modkit pileup $minimap2_primary_m4c/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam $ecoli_motif_focus_pileup/barcode${SLURM_ARRAY_TASK_ID}_m4c_gatcnnc.bed --motif GATCNNC 6 --ref $ecoli_ref

## targeting m5c
modkit pileup $minimap2_primary_m5c/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam $ecoli_motif_focus_pileup/barcode${SLURM_ARRAY_TASK_ID}_m5c_ccwgg.bed --motif CCWGG 1 --ref $ecoli_ref