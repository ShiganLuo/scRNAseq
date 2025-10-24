#!/bin/bash
function checkMd5() {
    md5Path=$1
    md5Dir=$(dirname "${md5Path}")
    md5File=$(basename "${md5Path}")
    cd ${md5Dir}
    md5sum -c ${md5File}
}
export -f checkMd5

dir=/disk5/luosg/scRNAseq/data

find "$dir" -type f -name "*md5.txt" | while read -r file; do
    echo "begin test ${file}"
    checkMd5 ${file}
done
