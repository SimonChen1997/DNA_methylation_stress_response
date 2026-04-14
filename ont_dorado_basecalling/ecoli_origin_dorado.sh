#!/bin/bash -l
#SBATCH --job-name="ecoli_origin_dorado"
#SBATCH --partition=gpu_cuda
#SBATCH --gres=gpu:h100:1
#SBATCH -o ecoli_origin_dorado.o
#SBATCH -e ecoli_origin_dorado.e

#########################################################
module use /sw/auto/rocky8c/epyc3_h100/modules/all/
module load cuda/12.6.0

#########################################################
### basecalling

$dorado basecaller --no-trim --recursive $model $input_scratch_pod5 --modified-bases 6mA 4mC_5mC \
    --kit-name SQK-NBD114-96 \
    > $output_scratch_bam/ecoli_origin.bam

$dorado demux --output-dir $output_scratch_demultiplex \
    --kit-name SQK-NBD114-96 $output_scratch_bam/ecoli_origin.bam