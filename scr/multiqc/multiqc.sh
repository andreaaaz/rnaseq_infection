#!/bin/bash

# 1. Asegurar que el sistema encuentre comandos instalados localmente en el HOME
export PATH=$HOME/.local/bin:$PATH

# 2. Instalar o verificar MultiQC localmente en tu cuenta
echo "=== Verificando / Instalando MultiQC ==="
pip install --user multiqc

# 3. Comprobar la versión instalada para el registro (log)
echo "=== Versión de MultiQC ==="
multiqc --version

# 4. Definir variables de ruta 
INPUT_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/exp/fastqc_results"
OUTPUT_DIR="/mnt/data/transcriptomica/iperez/proyectofinal/exp/multiqc_report"

# Crear la carpeta de destino si no existe
mkdir -p $OUTPUT_DIR

# 5. Ejecutar MultiQC
echo "=== Iniciando recopilación de reportes con MultiQC ==="
multiqc $INPUT_DIR -o $OUTPUT_DIR

echo "=== MultiQC finalizado con éxito ==="
