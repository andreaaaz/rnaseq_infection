#!/bin/bash
#SBATCH --job-name=fastqc_trimmed
#SBATCH --output=fastqc_trimmed_%j.log
#SBATCH --error=fastqc_trimmed_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4          # Usamos 4 cores para procesar en paralelo
#SBATCH --mem=16G                  
#SBATCH --time=02:00:00            # Suele ser más rápido con archivos trimmed

# 1. Cargar FastQC
module load fastqc/0.11.3

# 2. Definir rutas
# Apuntamos a la carpeta de salida de Cutadapt
INPUT_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/exp/trimmed_fastqs"
# Creamos una carpeta nueva para no mezclar estos resultados con los del fastqc raw
OUTPUT_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/exp/fastqc_trimmed_results"

mkdir -p $OUTPUT_DIR

# 3. Guardar todos los archivos recortados en una variable
FASTQ_FILES=$(ls $INPUT_DIR/*.fastq.gz)

echo "Archivos trimmed encontrados para procesar:"
echo "$FASTQ_FILES"
echo "-----------------------------------"

# 4. Correr FastQC optimizado con 4 hilos
fastqc -t 4 -o $OUTPUT_DIR $FASTQ_FILES

echo "¡FastQC de secuencias trimmed finalizado con éxito!"
