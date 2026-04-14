#!/bin/bash -l
#SBATCH --job-name="ecoli_origin_modkit_pileup"
#SBATCH --partition=general
#SBATCH --array=25-60
#SBATCH -o ecoli_origin_modkit_pileup.o
#SBATCH -e ecoli_origin_modkit_pileup.e

#########################################################
module load anaconda3
module load samtools

#########################################################
### extract fastq with methylation tag
samtools fastq -T ML,MM,MN $output_scratch_demultiplex/barcode${SLURM_ARRAY_TASK_ID}.bam > $fastq_sup/barcode${SLURM_ARRAY_TASK_ID}.fastq

#########################################################
### first round mapping to e.coli reference genome
source activate minimap2

minimap2 -y -ax lr:hq $ref $fastq_sup/barcode${SLURM_ARRAY_TASK_ID}.fastq -t 5 -o $minimap2/barcode${SLURM_ARRAY_TASK_ID}.bam

### index the sam file and convert to bam file
samtools view -bhS $minimap2/barcode${SLURM_ARRAY_TASK_ID}.bam | samtools sort -T $minimap2/barcode${SLURM_ARRAY_TASK_ID}.sorted -o $minimap2/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam
samtools index $minimap2/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

#########################################################
### extract only the primary mapped reads
samtools view -F 0x900 -F 4 -bhS $minimap2/barcode${SLURM_ARRAY_TASK_ID}.bam > $minimap2_primary/barcode${SLURM_ARRAY_TASK_ID}.bam

#########################################################
### extract fastq with methylation tag from primary bam
samtools fastq -T ML,MM,MN $minimap2_primary/barcode${SLURM_ARRAY_TASK_ID}.bam > $fastq_primary/barcode${SLURM_ARRAY_TASK_ID}.fastq

echo "${fastq_primary}/barcode${SLURM_ARRAY_TASK_ID}.fastq is done"

#########################################################
### remove short reads
source activate chopper

chopper -q 10 -l 250 -i $fastq_primary/barcode${SLURM_ARRAY_TASK_ID}.fastq > $chopper_dir/barcode${SLURM_ARRAY_TASK_ID}.fastq

echo "${chopper_dir}/barcode${SLURM_ARRAY_TASK_ID}.fastq"

#########################################################
### subsample fastq to 100x
source activate rasusa

rasusa reads -s 1727 --bases 1100MB $chopper_dir/barcode${SLURM_ARRAY_TASK_ID}.fastq -o $fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.fastq

#########################################################
### map to e.coli reference genome
source activate minimap2

minimap2 -y -ax lr:hq $ref $fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.fastq -t 5 -o $minimap2_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.bam

### index the sam file and convert to bam file
samtools view -bhS $minimap2_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.bam | \
samtools sort -T $minimap2_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted -o $minimap2_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam
samtools index $minimap2_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

#########################################################
### extract only the primary mapped reads
samtools view -F 0x900 -F 4 -bhS $minimap2_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.bam > $minimap2_primary_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.bam

#########################################################
### sort and index bam file
samtools view -bhS $minimap2_primary_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.bam | \
samtools sort -T $minimap2_primary_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted -o $minimap2_primary_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

samtools index $minimap2_primary_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

#########################################################
### adjust the modification information in bam file
source activate ont-modkit-0.5.0

## targeting m6a
modkit adjust-mods $minimap2_primary_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam stdout --ignore m | \
modkit adjust-mods stdin $minimap2_primary_m6a_200x/barcode${SLURM_ARRAY_TASK_ID}.bam --ignore 21839

## targeting m4c
modkit adjust-mods $minimap2_primary_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam stdout --ignore m | \
modkit adjust-mods stdin $minimap2_primary_m4c_200x/barcode${SLURM_ARRAY_TASK_ID}.bam --ignore a

## targeting m5c
modkit adjust-mods $minimap2_primary_fastq_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam stdout --ignore 21839 | \
modkit adjust-mods stdin $minimap2_primary_m5c_200x/barcode${SLURM_ARRAY_TASK_ID}.bam --ignore a

#########################################################
### sort and index bam file
## targeting m6a
samtools view -bhS $minimap2_primary_m6a_200x/barcode${SLURM_ARRAY_TASK_ID}.bam | \
samtools sort -T $minimap2_primary_m6a_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted -o $minimap2_primary_m6a_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

samtools index $minimap2_primary_m6a_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

## targeting m4c
samtools view -bhS $minimap2_primary_m4c_200x/barcode${SLURM_ARRAY_TASK_ID}.bam | \
samtools sort -T $minimap2_primary_m4c_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted -o $minimap2_primary_m4c_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

samtools index $minimap2_primary_m4c_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

## targeting m5c
samtools view -bhS $minimap2_primary_m5c_200x/barcode${SLURM_ARRAY_TASK_ID}.bam | \
samtools sort -T $minimap2_primary_m5c_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted -o $minimap2_primary_m5c_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

samtools index $minimap2_primary_m5c_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam

#########################################################
### get modification at base level
source activate ont-modkit-0.5.0

## targeting m6a
modkit pileup $minimap2_primary_m6a_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam $modkit_pileup_primary_m6a/barcode${SLURM_ARRAY_TASK_ID}.bed --motif A 0 --ref $ref

## targeting m4C
modkit pileup $minimap2_primary_m4c_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam $modkit_pileup_primary_m4c/barcode${SLURM_ARRAY_TASK_ID}.bed --motif C 0 --ref $ref

## targeting m5c
modkit pileup $minimap2_primary_m5c_200x/barcode${SLURM_ARRAY_TASK_ID}.sorted.bam $modkit_pileup_primary_m5c/barcode${SLURM_ARRAY_TASK_ID}.bed --motif C 0 --ref $ref
