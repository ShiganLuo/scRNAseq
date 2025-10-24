#!/bin/bash
Intestine=('CKO-chang-10XSC3' 'E2chang-10XSC3' 'TRA-chang-10XSC3' 'wt-chang-10XSC3')
Lung=('CKO-fei-10XSC3' 'E2fei-10XSC3' 'TRA-fei-10XSC3' 'wt-fei-10XSC3')
Muscle=('CKOF-10XSC3' 'E2jirou-10XSC3' 'TraF95-10XSC3' 'WT-10XSC3')
basedir=/ChIP_seq_2/StemCells/scRNAseq/output
purposedir=/disk5/luosg/scRNAseq/output/bam
for i in ${Intestine[@]};do
    # file=${basedir}/Intestine/${i}/outs/possorted_genome_bam.bam
    file=${basedir}/Intestine/${i}/outs/possorted_genome_bam.bam.bai
    mkdir -p ${purposedir}/Intestine/
    echo "process: ${file}"
    cp ${file} ${purposedir}/Intestine
    # mv ${purposedir}/Intestine/possorted_genome_bam.bam ${purposedir}/Intestine/${i}.bam
    mv ${purposedir}/Intestine/possorted_genome_bam.bam.bai ${purposedir}/Intestine/${i}.bam.bai
done
for i in ${Lung[@]};do
    # file=${basedir}/Lung/${i}/outs/possorted_genome_bam.bam
    file=${basedir}/Lung/${i}/outs/possorted_genome_bam.bam.bai
    mkdir -p ${purposedir}/Lung
    echo "process ${file}"
    cp ${file} ${purposedir}/Lung
    # mv ${purposedir}/Lung/possorted_genome_bam.bam ${purposedir}/Lung/${i}.bam
    mv ${purposedir}/Lung/possorted_genome_bam.bam.bai ${purposedir}/Lung/${i}.bam.bai
done
for i in ${Muscle[@]};do
    # file=/disk5/luosg/scRNAseq/output/algin/Muscle/${i}/outs/possorted_genome_bam.bam
    file=/disk5/luosg/scRNAseq/output/algin/Muscle/${i}/outs/possorted_genome_bam.bam.bai
    mkdir -p ${purposedir}/Muscle
    echo "process ${file}"
    cp ${file} ${purposedir}/Muscle
    # mv ${purposedir}/Muscle/possorted_genome_bam.bam ${purposedir}/Muscle/${i}.bam
    mv ${purposedir}/Muscle/possorted_genome_bam.bam.bai ${purposedir}/Muscle/${i}.bam.bai
done
