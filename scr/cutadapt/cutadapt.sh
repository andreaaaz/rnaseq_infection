#!/bin/bash
#SBATCH --job-name=cutadapt_trim
#SBATCH --output=cutadapt_%j.log
#SBATCH --error=cutadapt_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4          # 4 cores para acelerar la compresión y el recorte
#SBATCH --mem=12G                  
#SBATCH --time=04:00:00            

# 1. INSTALACIÓN LOCAL: Instalar Cutadapt en tu espacio de usuario
echo "Verificando/Instalando Cutadapt localmente..."
pip install --user cutadapt

# 2. Asegurar que el sistema encuentre el ejecutable de Cutadapt en tu home
export PATH=$HOME/.local/bin:$PATH

# Verificar que ya responda el comando
echo "Versión de Cutadapt instalada:"
cutadapt --version

# 3. Definir rutas del proyecto
BASE_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/data/fastqs"
OUTPUT_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/exp/trimmed_fastqs"

mkdir -p $OUTPUT_DIR

# 4. Secuencia del adaptador 
ADAPTER="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

echo "Iniciando el recorte de adaptadores..."

# 5. Iterar sobre cada carpeta de muestra
for dir in $BASE_DIR/SRR*/; do
    sample=$(basename "$dir")
    
    echo "--------------------------------------------"
    echo "Procesando muestra: $sample"
    
    INPUT_FILE="${dir}/${sample}_1.fastq"
    OUTPUT_FILE="${OUTPUT_DIR}/${sample}_trimmed.fastq.gz"
    
    # Ejecutar Cutadapt usando los 4 cores asignados
    cutadapt -a $ADAPTER -j 4 -m 20 -q 20 -o $OUTPUT_FILE $INPUT_FILE

done

echo "¡Procesamiento con Cutadapt finalizado con éxito!"
