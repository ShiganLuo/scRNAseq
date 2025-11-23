#!/bin/bash
#index
STAR=/opt/STAR-2.7.11b/bin/Linux_x86_64/STAR
function star_Index() {
    output_index=$1
    fasta=$2
    gtf=$3
    ${STAR} --runMode genomeGenerate \
    --runThreadN 50 \
    --genomeDir ${output_index} \
    --genomeFastaFiles ${fasta} \
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang 99
}
export -f star_Index
fasta=/home/lsg/Data/glioblastoma/data/ref/refdata-gex-GRCh38-2024-A/fasta/genome.fa
output_index=/home/lsg/Data/glioblastoma/output/index/star
gtf=/home/lsg/Data/glioblastoma/data/star/genes.gtf

#sjdbOverhang: max(ReadLength) -1;most comman sense 99 is proper
#comparsion

function star_align() {
    cat ${fq} | parallel -j 3 --colsep ' ' \
        ${STAR} --outSAMtype BAM SortedByCoordinate \
        --runThreadN 30 \
        --genomeDir ${output_index} \
        --readFilesIn {1} {2} \
        --outFileNamePrefix ${outdir}/{3} \
        --outSAMattributes NH HI AS nM CR CY UR UY \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 100 \
        --winAnchorMultimapNmax 100 \
        --outMultimapperOrder Random \
        --runRNGseed 777 \
        --outSAMmultNmax 1 \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist None \
        --soloBarcodeReadLength 98
}
fq=/home/lsg/Data/glioblastoma/output/index/fq.txt
outdir=/home/lsg/Data/glioblastoma/output/star


