#!/bin/bash
function renameFastq() {
    search_dir=$1
    # 查找所有 _002.fastq.gz 结尾的文件并重命名为 _001.fastq.gz
    find "$search_dir" -type f -name '*_002.fastq.gz' | while read filepath; do
        newpath="${filepath/_002.fastq.gz/_001.fastq.gz}"
        echo "Renaming:"
        echo "  From: $filepath"
        echo "  To:   $newpath"
        mv "$filepath" "$newpath"
    done
}
export -f renameFastq
