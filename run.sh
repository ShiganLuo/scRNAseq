#!/bin/bash

###################################  rename
source /disk5/luosg/scRNAseq/workflow/scRNAseq/scripts/check/rename.sh
# renameFastq /disk5/luosg/scRNAseq/data/fq

###################################  cellranger

source /disk5/luosg/scRNAseq/workflow/scRNAseq/scripts/1.o2.cellranger.sh

fqdir=/disk5/luosg/scRNAseq/data/fq
ref=/ChIP_seq_2/Data/index/Mus_musculus/CellRanger/GRCm39/refdata-gex-GRCm39-2024-A

Intestine=('CKO-chang-10XSC3' 'E2chang-10XSC3' 'TRA-chang-10XSC3' 'wt-chang-10XSC3')
Lung=('CKO-fei-10XSC3' 'E2fei-10XSC3' 'TRA-fei-10XSC3' 'wt-fei-10XSC3')
Muscle=('CKOF-10XSC3' 'E2jirou-10XSC3' 'TraF95-10XSC3' 'WT-10XSC3')

outdirPrefix=/disk5/luosg/scRNAseq/output/algin

# for sample in ${Intestine[@]:1:3};do
#     # echo ${sample}
#     fqdir1=${fqdir}/Intestine
#     mkdir -p ${outdirPrefix}/Intestine/
#     outdir=${outdirPrefix}/Intestine/${sample}
#     time cellrangerAlign ${fqdir1} ${sample} ${outdir} ${ref}
# done 

# for sample in ${Lung[@]:1:3};do
#     # echo ${sample}
#     fqdir1=${fqdir}/Lung
#     mkdir -p ${outdirPrefix}/Lung/
#     outdir=${outdirPrefix}/Lung/${sample}
#     time cellrangerAlign ${fqdir1} ${sample} ${outdir} ${ref}
# done 

# for sample in ${Muscle[@]};do
#     # echo ${sample}
#     fqdir1=${fqdir}/Muscle
#     mkdir -p ${outdirPrefix}/Muscle/
#     outdir=${outdirPrefix}/Muscle/${sample}
#     time cellrangerAlign ${fqdir1} ${sample} ${outdir} ${ref}
# done 

################################### bam filter

# Intestine=('CKO-chang-10XSC3' 'E2chang-10XSC3' 'TRA-chang-10XSC3' 'wt-chang-10XSC3')
# Lung=('CKO-fei-10XSC3' 'E2fei-10XSC3' 'TRA-fei-10XSC3' 'wt-fei-10XSC3')
# Muscle=('CKOF-10XSC3' 'E2jirou-10XSC3' 'TraF95-10XSC3' 'WT-10XSC3')

# basedir=/disk5/luosg/scRNAseq/output
# process_sample() {
#     group="$1"   
#     sample="$2"  

#     indir="${basedir}/bam/${group}"
#     outdir="${basedir}/bam_filter/${group}"
#     mkdir -p "${outdir}"
#     echo "processing:\n
#     samtools view ${indir}/${sample}.bam -h | awk '/^@/ || /CB:/' | samtools view -h -b > ${outdir}/${sample}_temp.bam"
#     samtools view ${indir}/${sample}.bam -h | awk '/^@/ || /CB:/' | samtools view -h -b > ${outdir}/${sample}_temp.bam
#     echo "processing:\n
#     samtools view ${outdir}/${sample}_temp.bam -h | awk '/^@/ || /UB:/' | samtools view -h -b > ${outdir}/${sample}_filter.bam"
#     samtools view ${outdir}/${sample}_temp.bam -h | awk '/^@/ || /UB:/' | samtools view -h -b > ${outdir}/${sample}_filter.bam
#     echo "processing:\n
#     rm -f ${outdir}/${sample}_temp.bam"
#     rm -f ${outdir}/${sample}_temp.bam
# }

# export -f process_sample
# export basedir

# groups=("Intestine" "Lung" "Muscle")

# task_list=()
# for group in "${groups[@]}"; do
#   declare -n samples="$group"      # samples 指向名为 $group 的数组
#   for sample in "${samples[@]}"; do
#     task_list+=("$group $sample")  # 组合为 "组名 样本名"
#   done
# done

# # 并行执行，每次运行 process_sample group sample
# printf '%s\n' "${task_list[@]}" | parallel -j 6 --colsep ' ' process_sample {1} {2}

################################### scTE
source /disk5/luosg/scRNAseq/workflow/scRNAseq/scripts/3.scTE.sh
# bam=$1
# out=$2
# index=$3
# thread=$4
index=/ChIP_seq_2/Data/index/Mus_musculus/scTE/mm39/mm39.exclusive.idx
thread=10
groups=("Intestine" "Lung" "Muscle")
basedir=/disk5/luosg/scRNAseq/output
task_list=()
for group in "${groups[@]}"; do
  declare -n samples="$group"      # samples 指向名为 $group 的数组
  for sample in "${samples[@]}"; do
    bam=${basedir}/bam_filter/${group}/${sample}_filter.bam
    outdir=${basedir}/scTE/${group}
    mkdir -p ${outdir}
    out=${basedir}/scTE/${group}/${sample}
    task_list+=("${bam} ${sample} ${outdir} ${index} ${thread}")  # 组合为 "组名 样本名"
  done
done
# echo ${task_list[@]}
# printf '%s\n' "${task_list[@]:1:11}" | parallel -j 1 --colsep ' ' scTEForCellranger {1} {2} {3} {4} {5}
printf '%s\n' "${task_list[@]:3:11}" | while read -r col1 col2 col3 col4 col5; do
    # echo "scTEForCellranger "$col1" "$col2" "$col3" "$col4" "$col5""
    scTEForCellranger "$col1" "$col2" "$col3" "$col4" "$col5"
done