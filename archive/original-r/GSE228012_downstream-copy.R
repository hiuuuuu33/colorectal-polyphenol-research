# ============================================================
# GSE228012 下游分析脚本
# 基于 GSE228012_RNA.R 的 limma 差异分析结果
# 参照 Mutigroups_GEOanalysisGSE54238.R 的分析模式
#
# 包含：
#   1. PCA 主成分分析
#   2. 火山图（3组对比）
#   3. 差异基因热图
#   4. GO 富集分析（BP / CC / MF）
#   5. KEGG 通路富集分析
# ============================================================
install.packages("AnnotationDbi")

# 清空环境
rm(list = ls())

# 加载包 ----------------------------------------------------------------
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(stringr)
library(ggrepel)

# PCA 相关
library(FactoMineR)
library(factoextra)

# 富集分析相关
library(clusterProfiler)
library(org.Hs.eg.db)

# 设置工作目录
setwd("/Users/hiugo/中药/数据分析")

# 加载上游分析结果 --------------------------------------------------------
load("GSE228012_limma_results.Rdata")

cat("\n")
cat("============================================================\n")
cat("  GSE228012 下游分析开始\n")
cat("============================================================\n\n")

# ============================================================
# 1. PCA 主成分分析
# ============================================================
cat(">>> 1. PCA 主成分分析 ...\n")

dat_pca <- as.data.frame(t(log2_fpkm))
dat.pca <- PCA(dat_pca, graph = FALSE)

pca_plot <- fviz_pca_ind(dat.pca,
                          geom.ind = "point",
                          col.ind = sample_info$group,
                          palette = "jco",
                          addEllipses = TRUE,
                          legend.title = "Groups",
                          title = "PCA Plot (GSE228012)",
                          axes = c(1, 2),
                         )


print(pca_plot)


ggsave("GSE228012_PCA.png", pca_plot, width = 8, height = 6)
cat("   PCA 图已保存至 GSE228012_PCA.png\n\n")

# ============================================================
# 2. 火山图 — 每个比较组分别绘制
# ============================================================
cat(">>> 2. 火山图 ...\n")

plot_volcano <- function(deg_full, contrast_name,
                          logFC_cutoff = 1,
                          pval_cutoff  = 0.05) {
  # deg_full 是 get_de_results 返回的 $full 数据框（含 gene、logFC、adj.P.Val 列）

  dat <- deg_full %>%
    mutate(
      neg_log10_pval = -log10(adj.P.Val),
      change = case_when(
        adj.P.Val < pval_cutoff & logFC >  logFC_cutoff ~ "Up",
        adj.P.Val < pval_cutoff & logFC < -logFC_cutoff ~ "Down",
        TRUE ~ "Stable"
      )
    )

  # 选取 top 基因标记
  top_up   <- dat %>% filter(change == "Up")   %>% arrange(desc(logFC)) %>% head(5)
  top_down <- dat %>% filter(change == "Down") %>% arrange(logFC)       %>% head(5)
  for_label <- bind_rows(top_up, top_down)

  p <- ggplot(dat, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(color = change), alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("Up"     = "#DC0000",
                                  "Down"   = "#3C5488",
                                  "Stable" = "grey70")) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
               lty = 4, col = "grey40", lwd = 0.6) +
    geom_hline(yintercept = -log10(pval_cutoff),
               lty = 4, col = "grey40", lwd = 0.6) +
    geom_label_repel(data = for_label,
                     aes(label = gene),
                     color = "black", size = 2.8,
                     max.overlaps = 20,
                     box.padding = 0.5) +
    labs(title    = paste0("Volcano: ", gsub("_", " ", contrast_name)),
         x        = "log2(Fold Change)",
         y        = "-log10(adjusted P-value)") +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank())

  filename <- paste0("GSE228012_volcano_", contrast_name, ".png")
  ggsave(filename, p, width = 9, height = 7)
  cat("   ", filename, "\n")
  return(p)
}



volcano_crc_normal     <- plot_volcano(crc_normal$full,     "CRC_vs_Normal")
volcano_adenoma_normal <- plot_volcano(adenoma_normal$full, "Adenoma_vs_Normal")
volcano_crc_adenoma    <- plot_volcano(crc_adenoma$full,    "CRC_vs_Adenoma")

cat("\n")


# ============================================================
# 3. 差异基因热图
# ============================================================
cat(">>> 3. 热图 ...\n")

# 取所有比较组显著差异基因的并集
all_sig_genes <- unique(c(
  rownames(crc_normal$significant),
  rownames(adenoma_normal$significant),
  rownames(crc_adenoma$significant)
))
cat("   所有比较组差异基因并集:", length(all_sig_genes), "\n")

# 若太多则截取 top 基因
if (length(all_sig_genes) > 2000) {
  top_crc_normal <- crc_normal$significant %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(abs(logFC))) %>%
    head(500)
  top_adenoma_normal <- adenoma_normal$significant %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(abs(logFC))) %>%
    head(500)
  top_crc_adenoma <- crc_adenoma$significant %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(abs(logFC))) %>%
    head(500)
  top_genes <- unique(c(top_crc_normal$gene,
                         top_adenoma_normal$gene,
                         top_crc_adenoma$gene))
  cat("   裁剪后基因数:", length(top_genes), "\n")
} else {
  top_genes <- all_sig_genes
}

heatmap_data <- log2_fpkm[top_genes, , drop = FALSE]
annotation_col <- data.frame(Group = sample_info$group,
                              row.names = colnames(heatmap_data))

# pheatmap 不是 ggplot 对象，用 png() + print() + dev.off()
png("GSE228012_heatmap.png", width = 10, height = 8, units = "in", res = 300)
heatmap_plot <- pheatmap(heatmap_data,
                          show_colnames   = FALSE,
                          show_rownames   = FALSE,
                          annotation_col  = annotation_col,
                          scale           = "row",
                          cluster_cols    = TRUE,
                          cluster_rows    = TRUE,
                          main            = paste0("DEGs Heatmap (", length(top_genes), " genes)\nGSE228012"))
dev.off()

cat("   GSE228012_heatmap.png 已保存\n\n")

# ============================================================
# 4. GO + KEGG 富集分析
# ============================================================
cat(">>> 4. GO / KEGG 富集分析 ...\n\n")

# 辅助函数：安全执行富集分析
safe_enrichGO <- function(gene_list, ont, ...) {
  tryCatch(
    enrichGO(gene          = gene_list,
             OrgDb         = org.Hs.eg.db,
             ont           = ont,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.05,
             readable      = TRUE,
             ...),
    error = function(e) {
      cat("    [!] GO", ont, "富集失败:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}

safe_enrichKEGG <- function(gene_list, ...) {
  tryCatch(
    enrichKEGG(gene         = gene_list,
               organism     = "hsa",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               ...),
    error = function(e) {
      cat("    [!] KEGG 富集失败:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}

# 富集分析主函数
run_enrichment <- function(sig_genes_df, contrast_name) {
  cat("  ---------- ", contrast_name, "----------\n")

  # 分离上下调基因
  up_symbols   <- rownames(sig_genes_df)[sig_genes_df$logFC > 0]
  down_symbols <- rownames(sig_genes_df)[sig_genes_df$logFC < 0]
  diff_symbols <- c(up_symbols, down_symbols)

  cat("    上调基因:", length(up_symbols), "\n")
  cat("    下调基因:", length(down_symbols), "\n")
  cat("    差异基因合计:", length(diff_symbols), "\n")

  if (length(diff_symbols) == 0) {
    cat("    [!] 无差异基因，跳过富集分析\n\n")
    return(NULL)
  }

  # Gene Symbol → ENTREZID
  s2e <- tryCatch(
    bitr(diff_symbols,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = org.Hs.eg.db),
    error = function(e) {
      cat("    [!] 基因 ID 转换失败:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(s2e) || nrow(s2e) == 0) {
    cat("    [!] 无可转换的 ENTREZID，跳过富集分析\n\n")
    return(NULL)
  }
  cat("    成功转换 ENTREZID:", nrow(s2e), "个基因\n")

  # 同时构建带 logFC 的有命名向量（用于 cnetplot 等）
  deg_with_id <- sig_genes_df %>%
    mutate(symbol = rownames(.)) %>%
    inner_join(s2e, by = c("symbol" = "SYMBOL"))
  geneList <- deg_with_id$logFC
  names(geneList) <- deg_with_id$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)

  # 上下调 ENTREZID
  up_id   <- s2e$ENTREZID[s2e$SYMBOL %in% up_symbols]
  down_id <- s2e$ENTREZID[s2e$SYMBOL %in% down_symbols]

  go_results   <- list()
  kegg_results <- list()

  # ---- GO 富集 (BP, CC, MF) ----
  for (ont in c("BP", "CC", "MF")) {
    go_res <- safe_enrichGO(s2e$ENTREZID, ont)
    if (!is.null(go_res) && nrow(go_res@result) > 0) {
      go_results[[ont]] <- go_res
      cat("    GO", ont, "显著条目:", nrow(go_res@result), "\n")

      # 气泡图
      dot_p <- dotplot(go_res, showCategory = 15,
                       title = paste0("GO ", ont, ": ", gsub("_", " ", contrast_name)))
      filename_png <- paste0("GSE228012_GO_", ont, "_", contrast_name, ".png")
      ggsave(filename_png, dot_p, width = 10, height = 8)

      # 输出 CSV
      write.csv(go_res@result,
                paste0("GSE228012_GO_", ont, "_", contrast_name, ".csv"),
                row.names = FALSE)
    }
  }

  # ---- KEGG 富集 (上调) ----
  if (length(up_id) >= 5) {
    kk_up <- safe_enrichKEGG(up_id)
    if (!is.null(kk_up) && nrow(kk_up@result) > 0) {
      kegg_results$up <- kk_up
      cat("    KEGG (上调) 显著条目:", nrow(kk_up@result), "\n")

      dot_p <- dotplot(kk_up, showCategory = 15,
                       title = paste0("KEGG Up: ", gsub("_", " ", contrast_name)))
      ggsave(paste0("GSE228012_KEGG_up_", contrast_name, ".png"),
             dot_p, width = 10, height = 8)
      write.csv(kk_up@result,
                paste0("GSE228012_KEGG_up_", contrast_name, ".csv"),
                row.names = FALSE)
    }
  } else {
    cat("    KEGG (上调) 跳过: 基因数不足5个\n")
  }

  # ---- KEGG 富集 (下调) ----
  if (length(down_id) >= 5) {
    kk_down <- safe_enrichKEGG(down_id)
    if (!is.null(kk_down) && nrow(kk_down@result) > 0) {
      kegg_results$down <- kk_down
      cat("    KEGG (下调) 显著条目:", nrow(kk_down@result), "\n")

      dot_p <- dotplot(kk_down, showCategory = 15,
                       title = paste0("KEGG Down: ", gsub("_", " ", contrast_name)))
      ggsave(paste0("GSE228012_KEGG_down_", contrast_name, ".png"),
             dot_p, width = 10, height = 8)
      write.csv(kk_down@result,
                paste0("GSE228012_KEGG_down_", contrast_name, ".csv"),
                row.names = FALSE)
    }
  } else {
    cat("    KEGG (下调) 跳过: 基因数不足5个\n")
  }

  # 保存该对比的富集结果
  save(s2e, go_results, kegg_results,
       file = paste0("GSE228012_enrichment_", contrast_name, ".Rdata"))

  cat("\n")
  return(list(s2e = s2e, go = go_results, kegg = kegg_results))
}

# 对 3 个比较组循环执行富集分析
enrich_crc_normal     <- run_enrichment(crc_normal$significant,     "CRC_vs_Normal")
enrich_adenoma_normal <- run_enrichment(adenoma_normal$significant, "Adenoma_vs_Normal")
enrich_crc_adenoma    <- run_enrichment(crc_adenoma$significant,    "CRC_vs_Adenoma")

# ============================================================
# 5. 保存全部结果
# ============================================================
cat(">>> 5. 保存全部结果 ...\n")

save(log2_fpkm, sample_info, fit2,
     crc_normal, adenoma_normal, crc_adenoma,
     pca_plot, volcano_crc_normal, volcano_adenoma_normal, volcano_crc_adenoma,
     heatmap_plot,
     enrich_crc_normal, enrich_adenoma_normal, enrich_crc_adenoma,
     file = "GSE228012_all_results.Rdata")

# ============================================================
# 6. 输出文件清单
# ============================================================
cat("\n")
cat("============================================================\n")
cat("  GSE228012 下游分析完成！\n")
cat("============================================================\n\n")
cat("  [PCA]       GSE228012_PCA.png\n")
cat("  [Volcano]   GSE228012_volcano_CRC_vs_Normal.png\n")
cat("  [Volcano]   GSE228012_volcano_Adenoma_vs_Normal.png\n")
cat("  [Volcano]   GSE228012_volcano_CRC_vs_Adenoma.png\n")
cat("  [Heatmap]   GSE228012_heatmap.png\n")
cat("  [GO/KEGG]   各比较组的 GO_*.png / GO_*.csv\n")
cat("  [GO/KEGG]   各比较组的 KEGG_*.png / KEGG_*.csv\n")
cat("  [Rdata]     GSE228012_enrichment_*.Rdata\n")
cat("  [Rdata]     GSE228012_all_results.Rdata\n")
cat("\n")
