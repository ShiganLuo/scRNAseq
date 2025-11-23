#!/bin/bash
scTE_build=/home/luosg/tools/scTE/bin/scTE_build
# index
# scTE_build -g hg38 -o /home/lsg/Data/glioblastoma/output/index/hg38 # Human
#########################################star#####################################

function scTEForStar() {
    scTE=/home/luosg/tools/scTE/bin/scTE
    local bam=$1
    local out=$2
    local outdir=$3
    local index=$4
    local thread=$5
    echo "excuete command: cd ${outdir}"
    cd ${outdir}
    echo "excute command: $scTE -i ${bam} -o ${out} -x ${index} -p ${thread} --min_counts 1 --min_genes 1 -CB CR -UMI UR"
    $scTE -i ${bam} -o ${out} -x ${index} -p ${thread} --min_counts 1 --min_genes 1 -CB CR -UMI UR
}
export -f scTEForStar
#########################################or cellranger#####################################

function scTEForCellranger() {
    scTE=/home/luosg/tools/scTE/bin/scTE
    local bam=$1
    local out=$2
    local outdir=$3
    local index=$4
    local thread=$5
    echo "excuete command: cd ${outdir}"
    cd ${outdir}
    echo "excute command: $scTE -i ${bam} -o ${out} -x ${index} -p ${thread} --min_counts 1 --min_genes 1 -CB CB -UMI UB"
    $scTE -i ${bam} -o ${out} -x ${index} -p ${thread} --min_counts 1 --min_genes 1 -CB CB -UMI UB
}
export -f scTEForCellranger

