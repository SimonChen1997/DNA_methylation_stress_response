#!/bin/bash -l
#SBATCH --job-name="ecoli_origin_condition_dmr"
#SBATCH --partition=general
#SBATCH -o ecoli_origin_condition_dmr.o
#SBATCH -e ecoli_origin_condition_dmr.e

#########################################################
module load anaconda3
module load rstudio/2024.04.2-r4.2.1

#########################################################
###### perform differntially methylation analysis using a premade R scripts
condition_array=("37" "42" "ACE" "PA" "RE" "TNA")

for phase in "ER" "SR"; do
    for i in ${!condition_array[@]}; do
        for j in $(seq $((i+1)) $((${#condition_array[@]}-1))); do
            if [[ "${condition_array[j]}" == "" || "${condition_array[j]}" == "${condition_array[i]}" ]]; then
                continue
            fi

            file_1=${condition_array[i]}_${phase}
            file_2=${condition_array[j]}_${phase}

            if [[ ${phase} == "ER" ]]; then
                growth_phase="exponential"
            elif [[ ${phase} == "SR" ]]; then
                growth_phase="stationary"
            fi

            echo "Now is comparing the file from ${phase}"
            echo "group a is ${condition_array[i]}"
            echo "group b is ${condition_array[j]}"

            mkdir $modkit_dmr_condition_summary/${condition_array[i]}_${condition_array[j]}_${phase}

### if investigating the origin and valid dataset correlation, $dmr_condition_parallel=dmr_condition_gene_pvalue.R
### if the differentially methylated sites within datasets, $dmr_condition_parallel=dmr_condition_parallel.R

            Rscript $dmr_condition_parallel \
                -p $modkit_dmr_condition \
                -d ecoli_m6a_${condition_array[i]}_${condition_array[j]}_${phase}_dmr_clean.tsv,ecoli_m4c_${condition_array[i]}_${condition_array[j]}_${phase}_dmr_clean.tsv,ecoli_m5c_${condition_array[i]}_${condition_array[j]}_${phase}_dmr_clean.tsv \
                -m $subset_condition_phase/${condition_array[i]}_${condition_array[j]}_m6a_${phase}_pileup.bed,$subset_condition_phase/${condition_array[i]}_${condition_array[j]}_m4c_${phase}_pileup.bed,$subset_condition_phase/${condition_array[i]}_${condition_array[j]}_m5c_${phase}_pileup.bed \
                -a ${condition_array[i]} \
                -b ${condition_array[j]} \
                -g $growth_phase \
                -o $modkit_dmr_condition_summary/${condition_array[i]}_${condition_array[j]}_${phase}

        done
    done
done
