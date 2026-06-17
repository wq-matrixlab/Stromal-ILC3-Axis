#!/bin/bash

GENOME_DIR="/data/QW/QW/scRNA_ref/Human_sc_refdata/STAR_ref/98"
WHITELIST="/data/QW/QW/scRNA_ref/Human_sc_refdata/STAR_ref/737K-august-2016.txt"
OUTPUT_BASE="/data/QW/QW/Public_Databases/GEO/Fibroblast_IBD_IL1/scVelo/"

SAMPLES=$(ls *R1_001.fastq.gz | sed 's/_L00.*//' | sort | uniq)

for SAMPLE in $SAMPLES; do
    mkdir -p ${OUTPUT_BASE}/${SAMPLE}/

    R1_FILES=$(ls ${SAMPLE}_L*_R1_001.fastq.gz | tr '\n' ',' | sed 's/,$//')
    R2_FILES=$(ls ${SAMPLE}_L*_R2_001.fastq.gz | tr '\n' ',' | sed 's/,$//')

    STAR \
        --runThreadN 16 \
        --genomeDir $GENOME_DIR \
        --readFilesIn $R2_FILES $R1_FILES \
        --readFilesCommand zcat \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist $WHITELIST \
        --soloCBstart 1 --soloCBlen 16 \
        --soloUMIstart 17 --soloUMIlen 10 \
        --soloFeatures Gene Velocyto \
        --soloBarcodeReadLength 0 \
        --outFileNamePrefix ${OUTPUT_BASE}/${SAMPLE}/ \
        --outSAMtype None
done