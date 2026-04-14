#!/bin/bash -l
#SBATCH --job-name="ecoli_origin_phase_dmr"
#SBATCH --partition=general
#SBATCH -o ecoli_origin_phase_dmr.o
#SBATCH -e ecoli_origin_phase_dmr.e

#########################################################
module load anaconda3
module load rstudio/2024.04.2-r4.2.1

#########################################################
###### perform differntially methylation analysis using a premade R scripts
for condition in "ACE" "TNA"; do

    echo "Now is comparing the phase file from ${condition}"
    mkdir $modkit_dmr_phase_gene/exp_sta_${condition}

    Rscript $dmr_phase_parallel \
        -p $modkit_dmr_phase \
        -d ecoli_m6a_exp_sta_${condition}_dmr_clean.tsv,ecoli_m4c_exp_sta_${condition}_dmr_clean.tsv,ecoli_m5c_exp_sta_${condition}_dmr_clean.tsv \
        -m $subset_phase/m6a_${condition}_pileup.bed,$subset_phase/m4c_${condition}_pileup.bed,$subset_phase/m5c_${condition}_pileup.bed \
        -a "exp" \
        -b "sta" \
        -c ${condition} \
        -s "origin" \
        -o $modkit_dmr_phase_gene/exp_sta_${condition}

done