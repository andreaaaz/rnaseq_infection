#!/bin/bash
#SBATCH --job-name=rnaseq_nf
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/rnaseq_nf_%j.log
#SBATCH --error=logs/rnaseq_nf_%j.err

module load nextflow

nextflow run rnaseq_infection.nf \
    --fastq "/mnt/data/transcriptomica/iperez/proyectofinal_nf/data/fastqs" \
    --bacteria "/mnt/data/transcriptomica/iperez/proyectofinal_nf/data/genomes/Escherichia_coli_K-12.fasta" \
    --phage "/mnt/data/transcriptomica/iperez/proyectofinal_nf/data/genomes/Enterobacteria_phage_T4.fasta" \
    --bacteria_annotation "/mnt/data/transcriptomica/iperez/proyectofinal_nf/data/genomes/Escherichia_coli_K-12.gff" \
    --phage_annotation "/mnt/data/transcriptomica/iperez/proyectofinal_nf/data/genomes/Enterobacteria_phage_T4.gff" \
    --outdir "./exp"

