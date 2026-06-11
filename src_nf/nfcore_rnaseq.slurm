#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=48:00:00

module load nextflow
module load singularity

nextflow run nf-core/rnaseq \
    -profile singularity \
    --input samplesheet.csv \
    --fasta genome.fasta \
    --gtf genes.gtf \
    --outdir results \
    --aligner hisat2 \
    -resume
