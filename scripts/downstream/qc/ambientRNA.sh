#!/bin/bash
cellbender remove-background \
    --input /disk5/luosg/scRNAseq/output/result/h5ad/Intestine_qc.h5ad \
    --output /disk5/luosg/scRNAseq/output/result/h5ad/Intestine_qc_removeAmbientRNA.h5ad \
    --cpu 5