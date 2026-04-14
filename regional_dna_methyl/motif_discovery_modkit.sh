#!/bin/bash -l
#SBATCH --job-name="motif_discovery_modkit"
#SBATCH --partition=general
#SBATCH --array=1-3
#SBATCH -o motif_discovery_modkit.o
#SBATCH -e motif_discovery_modkit.e

#########################################################
module load anaconda3

#########################################################
source activate ont-modkit-0.5.0

####### the exponential samples
for condition in "37" "42" "ACE" "PA" "RE" "TNA"; do

    echo "Now is comparing the phase file from ${condition} under exponential phase"

    modkit motif search -i $ecoli_clean_bed/${condition}_ER${SLURM_ARRAY_TASK_ID}_m6a_clean.bed -r $ecoli_ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 150 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $ecoli_motif_detection/${condition}_ER${SLURM_ARRAY_TASK_ID}_m6a_motifs.tsv --threads 20

    modkit motif search -i $ecoli_clean_bed/${condition}_ER${SLURM_ARRAY_TASK_ID}_m4c_clean.bed -r $ecoli_ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 150 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $ecoli_motif_detection/${condition}_ER${SLURM_ARRAY_TASK_ID}_m4c_motifs.tsv --threads 20

    modkit motif search -i $ecoli_clean_bed/${condition}_ER${SLURM_ARRAY_TASK_ID}_m5c_clean.bed -r $ecoli_ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 150 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $ecoli_motif_detection/${condition}_ER${SLURM_ARRAY_TASK_ID}_m5c_motifs.tsv --threads 20


done

####### the stationary samples
for condition in "37" "42" "ACE" "PA" "RE" "TNA"; do

    echo "Now is comparing the phase file from ${condition} under stationary phase"

    modkit motif search -i $ecoli_clean_bed/${condition}_SR${SLURM_ARRAY_TASK_ID}_m6a_clean.bed -r $ecoli_ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 150 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $ecoli_motif_detection/${condition}_SR${SLURM_ARRAY_TASK_ID}_m6a_motifs.tsv --threads 20

    modkit motif search -i $ecoli_clean_bed/${condition}_SR${SLURM_ARRAY_TASK_ID}_m4c_clean.bed -r $ecoli_ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 150 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $ecoli_motif_detection/${condition}_SR${SLURM_ARRAY_TASK_ID}_m4c_motifs.tsv --threads 20

    modkit motif search -i $ecoli_clean_bed/${condition}_SR${SLURM_ARRAY_TASK_ID}_m5c_clean.bed -r $ecoli_ref \
    --min-coverage 20 --low-thresh 0.01 --min-sites 150 --high-thresh 0.7 --min-frac-mod 0.7 --skip-search --min-log-odds 1.6 \
    -o $ecoli_motif_detection/${condition}_SR${SLURM_ARRAY_TASK_ID}_m5c_motifs.tsv --threads 20


done