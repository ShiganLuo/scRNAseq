#!/bin/bash
# wget -c -P /home/lsg/Data/glioblastoma/data/ref/ https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
# refMd5=$(md5sum /home/lsg/Data/glioblastoma/data/ref/refdata-gex-GRCh38-2024-A.tar.gz)
# echo -e "hg38基因组校验值\n${refMd5}"
function cellrangerAlign() {
    local cellranger=/opt/cellranger-9.0.1/bin/cellranger
    local fqdir=$1
    local sample=$2
    local outdir=$3
    local ref=$4
    $cellranger count \
        --fastqs=${fqdir} \
        --sample=${sample} \
        --transcriptome=${ref} \
        --id=${sample} \
        --output-dir=${outdir} \
        --create-bam=true \
        --localcores=35 \
        --nosecondary
}
export -f cellrangerAlign
# fqdir=/disk5/luosg/scRNAseq/data/fq
# sample='CKO-chang-10XSC3'
# outdir="/disk5/luosg/scRNAseq/output${sample}"
# ref=/ChIP_seq_2/Data/index/Mus_musculus/CellRanger/GRCm39/refdata-gex-GRCm39-2024-A
# cellrangerAlign ${fqdir} ${sample} ${outdir} ${ref}


