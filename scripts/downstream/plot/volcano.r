suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
volcano <- function(res, outjpeg, nlabel = 10, logFC = "logFC", Pvalue = "PValue", sort_by = "PValue") {
  res <- res %>%
    mutate(
      gene_status = case_when(
        logFC > 0.58 & -log10(PValue) > 1.3 ~ "Upregulated genes",
        logFC < -0.58 & -log10(PValue) > 1.3 ~ "Downregulated genes",
        TRUE ~ "Other genes"
      )
    )
  res <- res[!is.na(res$PValue), ]
  # Calculate the number of each gene_status
  status_counts <- res %>%
    group_by(gene_status) %>%
    summarise(count = n())

  # Create new legend labels containing quantities
  new_labels <- setNames(
    sapply(status_counts$gene_status, function(status) {
      count <- status_counts$count[status_counts$gene_status == status]
      paste0(status, " (n = ", count, ")")
    }),
    status_counts$gene_status
  )

  # get labels for the highest or lowest genes according to either PValue or
  # logFC
  if (sort_by == Pvalue) {
    up_genes <- res %>%
      filter(gene_status == "Upregulated genes") %>%
      arrange(PValue) %>%
      head(nlabel)
    down_genes <- res %>%
      filter(gene_status == "Downregulated genes") %>%
      arrange(PValue) %>%
      head(nlabel)
  } else if (sort_by == logFC) {
    up_genes <- head(arrange(
      res[res[, "gene_status"] == "Upregulated genes", ],
      desc(logFC)), nlabel)
    down_genes <- head(arrange(res[res[, "gene_status"] == "Upregulated genes", ], logFC), nlabel)
  } else {
    stop("Invalid sort_by argument. Choose either PValue or logFC.")
  }

  p <- ggplot(res, aes(logFC, -log10(PValue))) +
    geom_point(aes(size=-log10(PValue),color=-log10(PValue))) +
    scale_color_gradientn(values=seq(0,1,0.2),colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
    scale_size_continuous(range = c(1,3))+
    ggrepel::geom_text_repel(data = up_genes, aes(label = head(gene_name, nlabel)), color = "#F59494", size = 3) +
    ggrepel::geom_text_repel(data = down_genes, aes(label = head(gene_name, nlabel)), color = "#93ACF6", size = 3) +
    labs(x = expression(log[2]("FoldChange")), y = expression(-log[10]("adjusted p-value"))) +
    geom_vline(xintercept = 0.58, linetype = "dotted") +
    geom_vline(xintercept = -0.58, linetype = "dotted") +
    geom_hline(yintercept = 1.3, linetype = "dotted") +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  jpeg(width = 200, height = 150, units = "mm", res = 300, file = outjpeg)
  print(p)
  dev.off()
  print(paste0("volcano has bedd saved to", outjpeg))
}
# df = read.csv("/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_B_cell-combined_groupWT_chang_10XSC3_B_cell_DEG.tsv",sep="\t")
# df[['gene_name']] = rownames(df)
# print(head(df))
# rmsk = "/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz"
# df_rmsk <- read.table(rmsk, header = FALSE, sep = "\t")
# keep_vec = df_rmsk[[11]]
# df_filtered <- df[df$gene_name %in% keep_vec, ]
# print(head(df_rmsk[[11]]))
# volcano(df_filtered,"CKO_vs_WT.TE.jpeg",sort_by="PValue")
