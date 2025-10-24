from utils.TE import heatmap_DEG,heatmap_bulk, heatmap_sc
from utils.DEG_pseudo_bulk import Pre_AddAnnotations2RawCountsAnndata
import scanpy as sc
import logging
# outdir = "/disk5/luosg/scRNAseq/output/result/DEG/Intestine"
# outdir = "/disk5/luosg/scRNAseq/output/result/DEG/Lung"
outdir = "/disk5/luosg/scRNAseq/output/result/DEG/Muscle"
# adata = sc.read_h5ad("/disk5/luosg/scRNAseq/output/result/DEG/Intestine/h5ad/Intestine_bulk.h5ad")
# adata = sc.read_h5ad("/disk5/luosg/scRNAseq/output/result/DEG/Lung/h5ad/Lung_bulk.h5ad")
adata = sc.read_h5ad("/disk5/luosg/scRNAseq/output/result/DEG/Muscle/h5ad/Muscle_bulk.h5ad")
# adataRaw = sc.read_h5ad("/disk5/luosg/scRNAseq/output/result/h5ad/Intestine_qc.h5ad")
# adataAnn = sc.read_h5ad("/disk5/luosg/scRNAseq/output/result/h5ad/Intestine_annotate.h5ad")
# adata = Pre_AddAnnotations2RawCountsAnndata(adataRaw,adataAnn)
for cell in adata.obs["celltype"].unique():
    try:
        heatmap_bulk(
            adata,
            cell,
            tissue="jirou",
            out_image=f"{outdir}/plot/TE_heatmap_{cell}_bulk.png",
            rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz"
        )
        # heatmap_sc(
        #     adata,
        #     cell,
        #     f"{outdir}/plot/TE_heatmap_{cell}_sc.png",
        #     rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz"
        # )
    except ValueError as e:
        logging.warning(f"Skipping cell type '{cell}' due to an error: {e}")
        continue
    
logging.info("All available cell types have been processed.")
# heatmap_bulk(adata,"Enteroendocrine_cell",f"{outdir}/TE_heatmap_Enteroendocrine_cell_bulk.png",rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz")
# group_dict = {
#     "CKO_vs_WT":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_B_cell-combined_groupWT_chang_10XSC3_B_cell_DEG.tsv",
#     "E2_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupE2_chang_10XSC3_B_cell-combined_groupCKO_chang_10XSC3_B_cell_DEG.tsv",
#     "TRA_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupTRA_chang_10XSC3_B_cell-combined_groupCKO_chang_10XSC3_B_cell_DEG.tsv"
# }


# heatmap_DEG(group_dict,adata,"B_cell",out_image=f"{outdir}/TE_heatmap_B_cell_log2.png",rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz",method="log2")

# group_dict = {
#     "CKO_vs_WT":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_Enterocyte_cell-combined_groupWT_chang_10XSC3_Enterocyte_cell_DEG.tsv",
#     "E2_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupE2_chang_10XSC3_Enterocyte_cell-combined_groupCKO_chang_10XSC3_Enterocyte_cell_DEG.tsv",
#     "TRA_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupTRA_chang_10XSC3_Enterocyte_cell-combined_groupCKO_chang_10XSC3_Enterocyte_cell_DEG.tsv"
# }
# heatmap_DEG(group_dict,adata,"Enterocyte_cell",out_image=f"{outdir}/TE_heatmap_Enterocyte_cell.png",rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz")

# group_dict = {
#     "CKO_vs_WT":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_Enteroendocrine_cell-combined_groupWT_chang_10XSC3_Enteroendocrine_cell_DEG.tsv",
#     "E2_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupTRA_chang_10XSC3_Enteroendocrine_cell-combined_groupCKO_chang_10XSC3_Enteroendocrine_cell_DEG.tsv"
# }
# heatmap_DEG(group_dict,adata,"Enteroendocrine_cell",out_image=f"{outdir}/TE_heatmap_Enteroendocrine_cell.png",rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz")

# group_dict = {
#     "CKO_vs_WT":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_Intestinal_stem_cell-combined_groupWT_chang_10XSC3_Intestinal_stem_cell_DEG.tsv",
#     "E2_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupE2_chang_10XSC3_Intestinal_stem_cell-combined_groupCKO_chang_10XSC3_Intestinal_stem_cell_DEG.tsv",
#     "TRA_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupTRA_chang_10XSC3_Intestinal_stem_cell-combined_groupCKO_chang_10XSC3_Intestinal_stem_cell_DEG.tsv"
# }
# heatmap_DEG(group_dict,adata,"Intestinal_stem_cell",out_image=f"{outdir}/TE_heatmap_Intestinal_stem_cell.png",rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz")

# group_dict = {
#     "CKO_vs_WT":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_T_cytotoxic_cell-combined_groupWT_chang_10XSC3_T_cytotoxic_cell_DEG.tsv",
#     "E2_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupE2_chang_10XSC3_T_cytotoxic_cell-combined_groupCKO_chang_10XSC3_T_cytotoxic_cell_DEG.tsv",
#     "TRA_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupTRA_chang_10XSC3_T_cytotoxic_cell-combined_groupCKO_chang_10XSC3_T_cytotoxic_cell_DEG.tsv"
# }
# heatmap_DEG(group_dict,adata,"T_cytotoxic_cell",out_image=f"{outdir}/TE_heatmap_T_cytotoxic_cell.png",rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz")

# group_dict = {
#     "CKO_vs_WT":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_Tuft-combined_groupWT_chang_10XSC3_Tuft_DEG.tsv",
#     "TRA_vs_CKO":"/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupTRA_chang_10XSC3_Tuft-combined_groupCKO_chang_10XSC3_Tuft_DEG.tsv"
# }
# heatmap_DEG(group_dict,adata,"Tuft",out_image=f"{outdir}/TE_heatmap_Tuft.png",rmsk_file="/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz")


