#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%j.log
#SBATCH --error=fastqc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4          
#SBATCH --mem=16G                  
#SBATCH --time=03:00:00            

# 1. Cargar FastQC
module load fastqc/0.11.3

# 2. Definir rutas
BASE_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/data/fastqs"
OUTPUT_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/exp_"

mkdir -p $OUTPUT_DIR

# 3. Encontrar todos los archivos y guardarlos en una lista/variable
# 'find' buscará dentro de SRR21013136, SRR21013137, etc.
FASTQ_FILES=$(find $BASE_DIR -type f -name "*.fastq")

echo "Archivos encontrados para procesar:"
echo "$FASTQ_FILES"
echo "-----------------------------------"

# 4. Correr FastQC una sola vez pasándole toda la lista.
fastqc -t 4 -o $OUTPUT_DIR $FASTQ_FILES
