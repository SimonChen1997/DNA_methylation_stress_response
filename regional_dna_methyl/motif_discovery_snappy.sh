#!/bin/bash -l
#SBATCH --job-name="motif_discovery_snappy"
#SBATCH --partition=general
#SBATCH --array=1-3
#SBATCH -o motif_discovery_snappy.o
#SBATCH -e motif_discovery_snappy.e

#########################################################
module load anaconda3/2023.09-0

#########################################################
source activate snappy

for condition in "37" "42" "ACE" "PA" "RE" "TNA"; do

    ## targeting m6a
    snappy -mk_bed $modkit_pileup_primary_m6a/${condition}_ER${SLURM_ARRAY_TASK_ID}_m6a.bed -genome $sakai_ref -outdir $sakai_motif_snappy/${condition}_ER${SLURM_ARRAY_TASK_ID}_m6a

    ## targeting m4c
    snappy -mk_bed $modkit_pileup_primary_m4c/${condition}_ER${SLURM_ARRAY_TASK_ID}_m4c.bed -genome $sakai_ref -outdir $sakai_motif_snappy/${condition}_ER${SLURM_ARRAY_TASK_ID}_m4c

    ## targeting m5c
    snappy -mk_bed $modkit_pileup_primary_m5c/${condition}_ER${SLURM_ARRAY_TASK_ID}_m5c.bed -genome $sakai_ref -outdir $sakai_motif_snappy/${condition}_ER${SLURM_ARRAY_TASK_ID}_m5c

done

for condition in "37" "42" "ACE" "PA" "RE" "TNA"; do

    ## targeting m6a
    snappy -mk_bed $modkit_pileup_primary_m6a/${condition}_SR${SLURM_ARRAY_TASK_ID}_m6a.bed -genome $sakai_ref -outdir $sakai_motif_snappy/${condition}_SR${SLURM_ARRAY_TASK_ID}_m6a

    ## targeting m4c
    snappy -mk_bed $modkit_pileup_primary_m4c/${condition}_SR${SLURM_ARRAY_TASK_ID}_m4c.bed -genome $sakai_ref -outdir $sakai_motif_snappy/${condition}_SR${SLURM_ARRAY_TASK_ID}_m4c

    ## targeting m5c
    snappy -mk_bed $modkit_pileup_primary_m5c/${condition}_SR${SLURM_ARRAY_TASK_ID}_m5c.bed -genome $sakai_ref -outdir $sakai_motif_snappy/${condition}_SR${SLURM_ARRAY_TASK_ID}_m5c
    
done