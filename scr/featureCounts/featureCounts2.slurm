#!/bin/bash
#SBATCH --job-name=featureCounts2
#SBATCH --output=featureCounts_%j.log
#SBATCH --error=featureCounts_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# ------------------------------------------------------------------------------
# LOAD MODULES
# ------------------------------------------------------------------------------

module load subread/1.6.2

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------

PROJECT_DIR="/mnt/data/transcriptomica/iperez/proyectofinal"

BAM_DIR="${PROJECT_DIR}/exp/alignment_phage_bacteria_results/"
ANNOTATION="${PROJECT_DIR}/data/genomes/E_coli_and_T4_concat.gff"

OUT_DIR="${PROJECT_DIR}/exp/featureCounts"

mkdir -p $OUT_DIR

# ------------------------------------------------------------------------------
# FEATURECOUNTS
# ------------------------------------------------------------------------------

echo "===================================="
echo "Iniciando featureCounts"
echo "Fecha: $(date)"
echo "===================================="

featureCounts \
-T 24 \
-F GFF \
-t gene \
-g ID \
-s 1 \
-a ${ANNOTATION} \
-o ${OUT_DIR}/gene_counts.txt \
${BAM_DIR}/*.bam

echo "===================================="
echo "featureCounts terminado"
echo "Fecha: $(date)"
echo "===================================="
