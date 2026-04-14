#!/bin/bash -l
#SBATCH --job-name="ecoli_valid_modkit_dmr"
#SBATCH --partition=general
#SBATCH -o ecoli_valid_modkit_dmr.o
#SBATCH -e ecoli_valid_modkit_dmr.e
#SBATCH --account=a_qaafi_cas

#########################################################
module load anaconda3
module load samtools
module load htslib/1.15.1-gcc-11.3.0

#########################################################
### compress and index the pileup bed file
for file in $modkit_pileup_primary_motif/*.bed;do bgzip -k $file;done
for file in $modkit_pileup_primary_motif/*.bed.gz;do tabix -p bed $file;done

#########################################################
###### perform differntially methylation analysis
##### m6A for different conditions and phases
#### conditoins 
source activate ont-modkit-0.5.0
condition_array=("37" "ACEV" "TNAV")

for phase in "ER" "SR"; do
    for i in ${!condition_array[@]}; do
        for j in $(seq $((i+1)) $((${#condition_array[@]}-1))); do
            if [[ "${condition_array[j]}" == "" || "${condition_array[j]}" == "${condition_array[i]}" ]]; then
                continue
            fi

            file_1=${condition_array[i]}_${phase}
            file_2=${condition_array[j]}_${phase}

            echo "Now is comparing the file from ${phase}"
            echo "${file_1}1_m6a.bed"
            echo "${file_1}2_m6a.bed"
            echo "${file_1}3_m6a.bed"
            echo "and"
            echo "${file_2}1_m6a.bed"
            echo "${file_2}2_m6a.bed"
            echo "${file_2}3_m6a.bed"
            echo -e "\n \n"

            modkit dmr pair \
                -a $modkit_pileup_primary_motif/${file_1}1_m6a.bed.gz \
                -a $modkit_pileup_primary_motif/${file_1}2_m6a.bed.gz \
                -a $modkit_pileup_primary_motif/${file_1}3_m6a.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}1_m6a.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}2_m6a.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}3_m6a.bed.gz \
                --min-valid-coverage 20 \
                --ref $ref \
                --out-path $modkit_dmr_condition/ecoli_m6a_${condition_array[i]}_${condition_array[j]}_${phase}_dmr.tsv \
                --header \
                --base A \
                -t 4
        done
    done
done

#### phaes 
for condition in "37" "ACEV" "TNAV"; do
    exponential_file=${condition}_ER
    stationary_file=${condition}_SR

    echo "Now is comparing the file from ${condition}"
    echo "${exponential_file}1_m6a.bed"
    echo "${exponential_file}2_m6a.bed"
    echo "${exponential_file}3_m6a.bed"
    echo "and"
    echo "${stationary_file}1_m6a.bed"
    echo "${stationary_file}2_m6a.bed"
    echo "${stationary_file}3_m6a.bed"
    echo -e "\n \n"

    modkit dmr pair \
        -a $modkit_pileup_primary_motif/${exponential_file}1_m6a.bed.gz \
        -a $modkit_pileup_primary_motif/${exponential_file}2_m6a.bed.gz \
        -a $modkit_pileup_primary_motif/${exponential_file}3_m6a.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}1_m6a.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}2_m6a.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}3_m6a.bed.gz \
        --min-valid-coverage 20 \
        --ref $ref \
        --out-path $modkit_dmr_phase/ecoli_m6a_exp_sta_${condition}_dmr.tsv \
        --header \
        --base A \
        -t 4
done


##### m4C for different conditions and phases
#### conditoins 
condition_array=("37" "ACEV" "TNAV")

for phase in "ER" "SR"; do
    for i in ${!condition_array[@]}; do
        for j in $(seq $((i+1)) $((${#condition_array[@]}-1))); do
            if [[ "${condition_array[j]}" == "" || "${condition_array[j]}" == "${condition_array[i]}" ]]; then
                continue
            fi

            file_1=${condition_array[i]}_${phase}
            file_2=${condition_array[j]}_${phase}

            echo "Now is comparing the file from ${phase}"
            echo "${file_1}1_m4c.bed"
            echo "${file_1}2_m4c.bed"
            echo "${file_1}3_m4c.bed"
            echo "and"
            echo "${file_2}1_m4c.bed"
            echo "${file_2}2_m4c.bed"
            echo "${file_2}3_m4c.bed"
            echo -e "\n \n"

            modkit dmr pair \
                -a $modkit_pileup_primary_motif/${file_1}1_m4c.bed.gz \
                -a $modkit_pileup_primary_motif/${file_1}2_m4c.bed.gz \
                -a $modkit_pileup_primary_motif/${file_1}3_m4c.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}1_m4c.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}2_m4c.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}3_m4c.bed.gz \
                --min-valid-coverage 20 \
                --ref $ref \
                --out-path $modkit_dmr_condition/ecoli_m4c_${condition_array[i]}_${condition_array[j]}_${phase}_dmr.tsv \
                --header \
                --base C \
                -t 4
        done
    done
done

#### phaes 
for condition in "37" "ACEV" "TNAV"; do
    exponential_file=${condition}_ER
    stationary_file=${condition}_SR

    echo "Now is comparing the file from ${condition}"
    echo "${exponential_file}1_m4c.bed"
    echo "${exponential_file}2_m4c.bed"
    echo "${exponential_file}3_m4c.bed"
    echo "and"
    echo "${stationary_file}1_m4c.bed"
    echo "${stationary_file}2_m4c.bed"
    echo "${stationary_file}3_m4c.bed"
    echo -e "\n \n"

    modkit dmr pair \
        -a $modkit_pileup_primary_motif/${exponential_file}1_m4c.bed.gz \
        -a $modkit_pileup_primary_motif/${exponential_file}2_m4c.bed.gz \
        -a $modkit_pileup_primary_motif/${exponential_file}3_m4c.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}1_m4c.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}2_m4c.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}3_m4c.bed.gz \
        --min-valid-coverage 20 \
        --ref $ref \
        --out-path $modkit_dmr_phase/ecoli_m4c_exp_sta_${condition}_dmr.tsv \
        --header \
        --base C \
        -t 4
done

##### m5C for different conditions and phases
#### conditoins 
condition_array=("37" "ACEV" "TNAV")

for phase in "ER" "SR"; do
    for i in ${!condition_array[@]}; do
        for j in $(seq $((i+1)) $((${#condition_array[@]}-1))); do
            if [[ "${condition_array[j]}" == "" || "${condition_array[j]}" == "${condition_array[i]}" ]]; then
                continue
            fi

            file_1=${condition_array[i]}_${phase}
            file_2=${condition_array[j]}_${phase}

            echo "Now is comparing the file from ${phase}"
            echo "${file_1}1_m5c.bed"
            echo "${file_1}2_m5c.bed"
            echo "${file_1}3_m5c.bed"
            echo "and"
            echo "${file_2}1_m5c.bed"
            echo "${file_2}1_m5c.bed"
            echo "${file_2}3_m5c.bed"
            echo -e "\n \n"

            modkit dmr pair \
                -a $modkit_pileup_primary_motif/${file_1}1_m5c.bed.gz \
                -a $modkit_pileup_primary_motif/${file_1}2_m5c.bed.gz \
                -a $modkit_pileup_primary_motif/${file_1}3_m5c.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}1_m5c.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}2_m5c.bed.gz \
                -b $modkit_pileup_primary_motif/${file_2}3_m5c.bed.gz \
                --min-valid-coverage 20 \
                --ref $ref \
                --out-path $modkit_dmr_condition/ecoli_m5c_${condition_array[i]}_${condition_array[j]}_${phase}_dmr.tsv \
                --header \
                --base C \
                -t 4
        done
    done
done

#### phaes 
for condition in "37" "ACEV" "TNAV"; do
    exponential_file=${condition}_ER
    stationary_file=${condition}_SR

    echo "Now is comparing the file from ${condition}"
    echo "${exponential_file}1_m5c.bed"
    echo "${exponential_file}2_m5c.bed"
    echo "${exponential_file}3_m5c.bed"
    echo "and"
    echo "${stationary_file}1_m5c.bed"
    echo "${stationary_file}2_m5c.bed"
    echo "${stationary_file}3_m5c.bed"
    echo -e "\n \n"

    modkit dmr pair \
        -a $modkit_pileup_primary_motif/${exponential_file}1_m5c.bed.gz \
        -a $modkit_pileup_primary_motif/${exponential_file}2_m5c.bed.gz \
        -a $modkit_pileup_primary_motif/${exponential_file}3_m5c.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}1_m5c.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}2_m5c.bed.gz \
        -b $modkit_pileup_primary_motif/${stationary_file}3_m5c.bed.gz \
        --min-valid-coverage 20 \
        --ref $ref \
        --out-path $modkit_dmr_phase/ecoli_m5c_exp_sta_${condition}_dmr.tsv \
        --header \
        --base C \
        -t 4
done
