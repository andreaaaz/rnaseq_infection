#!/bin/bash
#SBATCH --job-name=HISAT2_phage
#SBATCH --output=hisat2_alignment_%j.out
#SBATCH --error=hisat2_alignment_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=danitzy2004@gmail.com
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --nodes=1

# loading modules
source /etc/profile.d/modules.sh
module load samtools/1.20
module load hisat2/2.2.1

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------

REF_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/data/genomes"
REF_FASTA="${REF_DIR}/E_coli_and_T4_concat.fasta"
INDEX_BASE="${REF_DIR}/hisat2_index/E_coli_and_T4"

THREADS=24

# ------------------------------------------------------------------------------
# INDEX REFERENCE GENOME
# ------------------------------------------------------------------------------
if [ ! -f "${INDEX_BASE}.1.ht2" ]; then
    echo "=== Creando índices de HISAT2 ==="
    mkdir -p "${REF_DIR}/hisat2_index"
    hisat2-build -p $THREADS $REF_FASTA $INDEX_BASE
else
    echo "=== Los índices ya existen. ==="
fi

# ------------------------------------------------------------------------------
# Single-End ALIGNMENT
# ------------------------------------------------------------------------------
echo "=== Iniciando mapeo con HISAT2 ==="

#ruta donde están FASTQs
FASTQ_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/exp/trimmed_fastqs"

for file in $(find "$FASTQ_DIR" -maxdepth 1 -name "*_trimmed.fastq.gz")
do
    #Extraer nombre del archivo
    filename=$(basename "$file")

    #Extraer el ID de la muestra
    base=$(basename "$filename" _trimmed.fastq.gz)

    echo "Procesando muestra: $base"

    #Mapeo con HISAT2
    # Mapeo, conversión y ordenamiento directo
    hisat2 -p $THREADS -x $INDEX_BASE -U "$file" 2> ${base}_hisat2.log | \
    samtools sort -@ $THREADS -m 2G -o ${base}.bam -

    echo "  Indexando BAM..."
    samtools index ${base}.bam

    echo "Muestra $base finalizada exitosamente."
    echo "------------------------------------------------"
done

echo "=== Alineamineto terminado ==="
exit;

