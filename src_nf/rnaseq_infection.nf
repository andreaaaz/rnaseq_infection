#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Paramateros

params.fastq = "./data/fastqs/*.fastq.gz"
params.bacteria = "./data/genomes/E_coli.fa"
params.phage = "./data/genomes/T4.fa"
params.phage_annotation ="./data/genomes/T4.gff"
params.bacteria_annotation ="./data/genomes/E_coli.gff"
params.outdir = "./exp"
params.adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

// Canales

Channel
    .fromPath(params.fastq)
    .map { file ->
        tuple(file.baseName,file)
    }
    .set { raw_fastqs }

bacteria_ch = Channel.fromPath(params.bacteria)
phage_ch = Channel.fromPath(params.phage)
bac_ann_ch = Channel.fromPath(params.bacteria_annotation)
phag_ann_ch = Channel.fromPath(params.phage_annotation)

// Quality control

process FASTQC_RAW {

    tag "$sample"
    publishDir "${params.outdir}/quality1", mode: 'copy'

    input:
    tuple val(sample), path(fastq)
    
    output:
    tuple val(sample), path(fastq)
    path("*_fastqc.html")
    path("*_fastqc.zip")

    script:
    """
    fastqc \
        -t ${task.cpus} \
        -o . \
        ${fastq}
    """
}


process CUTADAPT {

    tag "$sample"
    publishDir "${params.outdir}/TRIM_results", mode: 'copy'

    input:
    tuple val(sample), path(fastq)
    
    output:
    tuple(val(sample), path("${sample}_trimmed.fastq.gz"))

    script:
    """
    cutadapt \
        -a ${params.adapter} \
        -j ${task.cpus} \
        -m 20 \
        -q 20 \
        -o ${sample}_trimmed.fastq.gz \
        ${fastq}
    """
}

process FASTQC_TRIMMED {

    tag "$sample"
    publishDir"${params.project_dir}/exp/quality2", mode: 'copy'

    input:
    tuple val(sample), path(fastq)
    
    output:
    tuple val(sample), path(fastq)
    path("*_fastqc.html") 
    path("*_fastqc.zip")

    script:
    """
    fastqc \
        -t ${task.cpus} \
        -o . \
        ${fastq}
    """
}

// Alineamiento

process CONCAT {

    publishDir( "${params.outdir}/reference", mode: 'copy' )

    input:
    path(bacteria)
    path(phage)
    path(bac_ann)
    path(phag_ann)

    output:
    tuple path("E_coli_and_T4.fasta")
    path("E_coli_and_T4.gff")
    
    script:
    """
    cat ${bacteria} ${phage} > E_coli_and_T4.fasta
    cat ${bac_ann} ${phag_ann} > E_coli_and_T4.gff
    """
}

process HISAT2_BUILD {

    publishDir "${params.outdir}/hisat2_index", mode: 'copy'

    input:
    tuple path(fasta), path(gff)
    
    output:
    path("hisat2_index")

    script:
    """
    mkdir hisat2_index
    hisat2-build \
        -p ${task.cpus} \
        ${fasta} \
        hisat2_index/E_coli_and_T4
    """
}

process BOWTIE2_BUILD {

    publishDir "${params.outdir}/bowtie2_index", mode: 'copy'

    input:
    tuple path(fasta), path(gff)
    
    output:
    path("bowtie2_index")

    script:
    """
    mkdir bowtie2_index
    bowtie2-build \
        ${fasta} \
        bowtie2_index/E_coli_and_T4
    """
}

process HISAT2_ALIGN {

    tag "$sample"
    publishDir "${params.outdir}/hisat2_alignment", mode: 'copy'
    
    input:
    path(index_dir)
    tuple val(sample), path(fastq)
    
    output:
    tuple val(sample), 
          path("${sample}.bam"), 
          path("${sample}.bam.bai")

    script:
    """
    hisat2 \
        -p ${task.cpus} \
        -x ${index_dir}/E_coli_and_T4 \
        -U ${fastq} \
        2> ${sample}_hisat2.log | \
    samtools sort \
        -@ ${task.cpus} \
        -o ${sample}.bam -
    samtools index ${sample}.bam
    """
}


process BOWTIE2_ALIGN {

    tag "$sample"
    publishDir"${params.outdir}/bowtie2_alignment",mode: 'copy'

    input:
    path(index_dir)
    tuple val(sample), path(fastq)
    
    output:
    tuple val(sample),
          path("${sample}.bam"),
          path("${sample}.bam.bai")

    script:
    """
    bowtie2 \
        -p ${task.cpus} \
        -x ${index_dir}/E_coli_and_T4 \
        -U ${fastq} \
        2> ${sample}_bowtie2.log | \
    samtools sort \
        -@ ${task.cpus} \
        -o ${sample}.bam -
    samtools index ${sample}.bam
    """
}


// Features counts

process FEATURECOUNTS_HISAT2 {

    publishDir "${params.outdir}/featureCounts_hisat2", mode: 'copy'

    input:
    tuple path(fasta), path(gff)
    path(bams)
    
    output:
    path("gene_counts_hisat2.txt")

    script:
    def bam_string = bams.join(' ')
    """
    featureCounts \
        -T ${task.cpus} \
        -F GFF \
        -t gene \
        -g ID \
        -s 2 \
        -a ${gff} \
        -o gene_counts_hisat2.txt \
        ${bam_string}
    """
}

process FEATURECOUNTS_BOWTIE2 {

    publishDir "${params.outdir}/featureCounts_bowtie2", mode: 'copy'

    input:
    tuple path(fasta), path(gff)
    path(bams)
    
    output:
    path("gene_counts_bowtie2.txt")
    
    script:
    def bam_string = bams.join(' ')
    """
    featureCounts \
        -T ${task.cpus} \
        -F GFF \
        -t gene \
        -g ID \
        -s 2 \
        -a ${params.annotation} \
        -o gene_counts_bowtie2.txt \
        ${bam_string}
    """
}

process MULTIQC {

    publishDir"${params.outdir}/multiqc", mode: 'copy'

    input:
    path(files)
    
    output:
    path("multiqc_report.html")

    script:
    """
    multiqc .  -o .
    """
}

// Workflow

workflow {

    FASTQC_RAW(raw_fastqs)
    trimmed = CUTADAPT(raw_fastqs)
    FASTQC_TRIMMED(trimmed)
    trimmed.into {
        trimmed_hisat2
        trimmed_bowtie2
    }

    reference = CONCAT(bacteria_ch, phage_ch, bac_ann_ch, phag_ann_ch,)
    
    hisat2_index = HISAT2_BUILD(reference)
    hisat2_bams = HISAT2_ALIGN(hisat2_index,trimmed_hisat2)
    hisat2_bam_files = hisat2_bams.map {
      sample, bam, bai -> bam
    }
    
    FEATURECOUNTS_HISAT2(
      reference,
      hisat2_bam_files.collect()
    )

    bowtie2_index = BOWTIE2_BUILD(reference)
    bowtie2_bams = BOWTIE2_ALIGN(bowtie2_index,trimmed_bowtie2)
    bowtie2_bam_files = bowtie2_bams.map {
      sample, bam, bai -> bam
    }

    FEATURECOUNTS_BOWTIE2(
      reference,
      bowtie2_bam_files.collect()
    )
}
