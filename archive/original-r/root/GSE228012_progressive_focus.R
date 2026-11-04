# ============================================================
# GSE228012 全程上调/下调基因分析
#
# 核心思路：
#   从差异分析显著性出发，取在 Normal→Adenoma 和 Adenoma→CRC
#   两个阶段中方向一致且均显著的基因，定义为"全程上调/下调基因"
#
# 分析内容：
#   1. 全程上调/下调基因鉴定
#   2. 代表性折线图（标注基因名）
#   3. 全程基因表达热图
#   4. 与 CRC_vs_Normal 的 Venn 重叠分析
#   5. 全程基因 GO / KEGG 富集分析
#   6. 已知 CRC 驱动基因的跨阶段表达
# ============================================================

rm(list = ls())

# ---- 加载包 ----
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggrepel)

setwd("d:/Bolejiaoyu/maixiao")

# ---- 创建输出文件夹 ----
outdir <- "progressive_analysis"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
cat("输出目录:", outdir, "\n")

# ---- 加载上游结果 ----
load("GSE228012_limma_results.Rdata")

cat("\n============================================================\n")
cat("  GSE228012 全程上调/下调基因分析\n")
cat("============================================================\n\n")

# ============================================================
# 1. 鉴定全程上调/下调基因
# ============================================================
cat(">>> 1. 鉴定全程上调/下调基因 ...\n")

# 从两个过渡阶段的显著基因中提取
na_up   <- rownames(adenoma_normal$significant)[adenoma_normal$significant$logFC > 0]
na_down <- rownames(adenoma_normal$significant)[adenoma_normal$significant$logFC < 0]
ac_up   <- rownames(crc_adenoma$significant)[crc_adenoma$significant$logFC > 0]
ac_down <- rownames(crc_adenoma$significant)[crc_adenoma$significant$logFC < 0]

# 全程上调：在 N→A 上调 且 A→CRC 上调
prog_up <- intersect(na_up, ac_up)
# 全程下调：在 N→A 下调 且 A→CRC 下调
prog_down <- intersect(na_down, ac_down)

cat("  N→A 上调基因:      ", length(na_up), "\n")
cat("  A→CRC 上调基因:    ", length(ac_up), "\n")
cat("  全程上调基因 (交集):", length(prog_up), "\n\n")

cat("  N→A 下调基因:      ", length(na_down), "\n")
cat("  A→CRC 下调基因:    ", length(ac_down), "\n")
cat("  全程下调基因 (交集):", length(prog_down), "\n\n")

# 构建全程基因综合表（合并两个阶段的 logFC 和显著性）
build_prog_table <- function(gene_list, direction) {
  if (length(gene_list) == 0) return(NULL)

  na_res <- adenoma_normal$full
  ac_res <- crc_adenoma$full
  nc_res <- crc_normal$full

  m_na <- match(gene_list, na_res$gene)
  m_ac <- match(gene_list, ac_res$gene)
  m_nc <- match(gene_list, nc_res$gene)

  data.frame(
    gene        = gene_list,
    logFC_NA    = na_res$logFC[m_na],
    adjP_NA     = na_res$adj.P.Val[m_na],
    logFC_AC    = ac_res$logFC[m_ac],
    adjP_AC     = ac_res$adj.P.Val[m_ac],
    logFC_NC    = nc_res$logFC[m_nc],
    adjP_NC     = nc_res$adj.P.Val[m_nc],
    direction   = direction,
    stringsAsFactors = FALSE
  ) %>% arrange(adjP_NA)  # 按第一阶段的显著性排序
}

prog_up_table   <- build_prog_table(prog_up,   "Up-Up")
prog_down_table <- build_prog_table(prog_down, "Down-Down")

# 合并
prog_all <- bind_rows(prog_up_table, prog_down_table)
cat("  全程上调基因: ", nrow(prog_up_table), "\n")
cat("  全程下调基因: ", nrow(prog_down_table), "\n")
cat("  全程变化基因合计:", nrow(prog_all), "\n\n")

# 输出 CSV
write.csv(prog_all, file.path(outdir, "progressive_genes.csv"), row.names = FALSE)
cat("  progressive_genes.csv 已保存\n\n")

# ============================================================
# 2. 代表性基因拆线图（标注基因名）
# ============================================================
cat(">>> 2. 代表性基因拆线图 ...\n")

# 按 |logFC_NC| 排序，各取 top 15
top_up   <- prog_up_table   %>% arrange(desc(abs(logFC_NC))) %>% head(15)
top_down <- prog_down_table %>% arrange(desc(abs(logFC_NC))) %>% head(15)
top_all  <- bind_rows(top_up, top_down)

# 计算每个基因在三组中的平均表达
expr_long <- data.frame()
for (i in seq_len(nrow(top_all))) {
  g <- top_all$gene[i]
  if (!g %in% rownames(log2_fpkm)) next
  vals <- log2_fpkm[g, ]
  expr_long <- rbind(expr_long, data.frame(
    gene      = g,
    direction = top_all$direction[i],
    Normal    = mean(as.numeric(vals[sample_info$group == "Normal"]),  na.rm = TRUE),
    Adenoma   = mean(as.numeric(vals[sample_info$group == "Adenoma"]), na.rm = TRUE),
    CRC       = mean(as.numeric(vals[sample_info$group == "CRC"]),     na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}

# 转为长格式
expr_long <- expr_long %>%
  pivot_longer(cols = c(Normal, Adenoma, CRC),
               names_to = "Stage", values_to = "Expression") %>%
  mutate(Stage = factor(Stage, levels = c("Normal", "Adenoma", "CRC")))

# 分段绘制（上调一组，下调一组），确保标签不重叠
plot_lines <- function(data, dir_label, color_val, filename) {
  p <- ggplot(data, aes(x = Stage, y = Expression, group = gene)) +
    geom_line(aes(color = gene), linewidth = 1, alpha = 0.8) +
    geom_point(aes(color = gene), size = 2.5) +
    geom_label_repel(
      data = data %>% filter(Stage == "CRC"),
      aes(label = gene, color = gene),
      size = 3.2,
      direction = "y",
      nudge_x = 0.4,
      segment.size = 0.3,
      show.legend = FALSE,
      max.overlaps = 20,
      box.padding = 0.4
    ) +
    scale_color_manual(values = setNames(
      rep(color_val, length(unique(data$gene))),
      unique(data$gene)
    ), guide = "none") +
    labs(title = paste0("Progressive ", dir_label, " Genes (Top 15)"),
         subtitle = "Normal → Adenoma → CRC",
         x = "", y = "Mean log2(FPKM + 1)") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p)
  ggsave(file.path(outdir, filename), p, width = 11, height = 8)
  cat("  ", filename, " 已保存\n")
}

plot_lines(
  expr_long %>% filter(direction == "Up-Up"),
  "Up-regulated", "#DC0000",
  "progressive_up_lines.png"
)

plot_lines(
  expr_long %>% filter(direction == "Down-Down"),
  "Down-regulated", "#3C5488",
  "progressive_down_lines.png"
)

# ============================================================
# 3. 全程基因热图
# ============================================================
cat("\n>>> 3. 全程基因表达热图 ...\n")

prog_genes <- unique(c(prog_up, prog_down))
heat_mat <- log2_fpkm[prog_genes, , drop = FALSE]
heat_z   <- t(scale(t(heat_mat)))  # z-score

# 样本按分组排序
stage_ord <- order(sample_info$group)
heat_z    <- heat_z[, stage_ord]

# 行注释
row_anno <- data.frame(
  Direction = ifelse(prog_genes %in% prog_up, "Up-Up", "Down-Down"),
  row.names = prog_genes
)
# 列注释
col_anno <- data.frame(
  Stage = sample_info$group[stage_ord],
  row.names = colnames(heat_z)
)

ann_colors <- list(
  Direction = c("Up-Up" = "#DC0000", "Down-Down" = "#3C5488"),
  Stage     = c(Normal = "#00A087", Adenoma = "#E64B35", CRC = "#3C5488")
)

# 显示的基因标签：取 top 30 (各 15)
show_labels <- c(
  prog_up_table$gene[1:min(15, nrow(prog_up_table))],
  prog_down_table$gene[1:min(15, nrow(prog_down_table))]
)
rownames(heat_z)[!rownames(heat_z) %in% show_labels] <- ""

png(file.path(outdir, "progressive_heatmap.png"),
    width = 12, height = 10, units = "in", res = 300)
ph <- pheatmap(heat_z,
  show_colnames   = FALSE,
  show_rownames   = TRUE,
  fontsize_row    = 6,
  annotation_col  = col_anno,
  annotation_row  = row_anno,
  annotation_colors = ann_colors,
  scale           = "none",
  cluster_cols    = FALSE,
  cluster_rows    = TRUE,
  gaps_col        = cumsum(table(col_anno$Stage)),
  main            = paste0("Progressive Genes Heatmap (", nrow(heat_z), " genes)\nNormal → Adenoma → CRC"))
print(ph)
dev.off()
cat("  progressive_heatmap.png 已保存\n\n")

# ============================================================
# 4. 重叠分析（Venn：全程上调 vs CRC_vs_Normal 上调）
# ============================================================
cat(">>> 4. 重叠分析 ...\n")

# 定义各基因集合
nc_up   <- rownames(crc_normal$significant)[crc_normal$significant$logFC > 0]
nc_down <- rownames(crc_normal$significant)[crc_normal$significant$logFC < 0]

gene_sets <- list(
  Progressive_Up   = prog_up,
  Progressive_Down = prog_down,
  CRC_vs_Normal_Up = nc_up,
  CRC_vs_Normal_Down = nc_down
)

for (nm in names(gene_sets)) {
  cat(sprintf("  %-22s: %d\n", nm, length(gene_sets[[nm]])))
}

# ---- Venn 上调：Progressive_Up vs NA_Up vs AC_Up ----
cat("\n  绘制 Venn / 重叠图 ...\n")

# 上调用 Venn
png(file.path(outdir, "venn_progressive_up.png"),
    width = 8, height = 8, units = "in", res = 200)
up_list <- list(
  `N→A Up`  = na_up,
  `A→CRC Up`= ac_up,
  `N→CRC Up`= nc_up
)
up_counts <- vennCounts(
  cbind(
    as.integer(unique(c(na_up, ac_up, nc_up)) %in% na_up),
    as.integer(unique(c(na_up, ac_up, nc_up)) %in% ac_up),
    as.integer(unique(c(na_up, ac_up, nc_up)) %in% nc_up)
  )
)
vennDiagram(up_counts,
  main = "Up-regulated Gene Overlap\nAcross Stage Transitions",
  circle.col = c("#E64B35", "#F39B7F", "#DC0000"))
dev.off()
cat("  venn_progressive_up.png 已保存\n")

# 下调用 Venn
png(file.path(outdir, "venn_progressive_down.png"),
    width = 8, height = 8, units = "in", res = 200)
down_list <- list(
  `N→A Down`  = na_down,
  `A→CRC Down`= ac_down,
  `N→CRC Down`= nc_down
)
down_counts <- vennCounts(
  cbind(
    as.integer(unique(c(na_down, ac_down, nc_down)) %in% na_down),
    as.integer(unique(c(na_down, ac_down, nc_down)) %in% ac_down),
    as.integer(unique(c(na_down, ac_down, nc_down)) %in% nc_down)
  )
)
vennDiagram(down_counts,
  main = "Down-regulated Gene Overlap\nAcross Stage Transitions",
  circle.col = c("#4DBBD5", "#8491B4", "#3C5488"))
dev.off()
cat("  venn_progressive_down.png 已保存\n\n")

# ============================================================
# 5. GO / KEGG 富集分析
# ============================================================
cat(">>> 5. GO / KEGG 富集分析 ...\n\n")

# 延迟加载 Bioconductor 包，避免 dplyr 屏蔽
library(clusterProfiler)
library(org.Hs.eg.db)

# ID 转换
convert_id <- function(genes) {
  tryCatch({
    s2e <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
    return(s2e$ENTREZID)
  }, error = function(e) {
    cat("    [!] ID转换失败:", conditionMessage(e), "\n")
    return(character(0))
  })
}

# 安全富集
do_go <- function(entrez, ont_label = "BP") {
  if (length(entrez) < 5) return(NULL)
  tryCatch(
    enrichGO(gene = entrez, OrgDb = org.Hs.eg.db,
             ont = ont_label, pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.05,
             readable = TRUE),
    error = function(e) { cat("    [!] GO失败:", conditionMessage(e), "\n"); NULL }
  )
}

do_kegg <- function(entrez) {
  if (length(entrez) < 5) return(NULL)
  tryCatch(
    enrichKEGG(gene = entrez, organism = "hsa",
               pvalueCutoff = 0.05, qvalueCutoff = 0.05),
    error = function(e) { cat("    [!] KEGG失败:", conditionMessage(e), "\n"); NULL }
  )
}

# 对全程上调和全程下调分别做富集
enrich_list <- list(
  Progressive_Up   = prog_up,
  Progressive_Down = prog_down
)

for (grp in names(enrich_list)) {
  genes <- enrich_list[[grp]]
  cat("  ---", grp, "---\n")
  cat("    基因数:", length(genes), "\n")

  if (length(genes) < 5) { cat("    基因数不足，跳过\n\n"); next }

  entrez <- convert_id(genes)
  cat("    ENTREZID:", length(entrez), "\n")

  # GO BP
  go_bp <- do_go(entrez, "BP")
  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    cat("    GO BP:", nrow(go_bp@result), "terms\n")
    tryCatch({
      p_bp <- dotplot(go_bp, showCategory = min(15, nrow(go_bp@result)),
                      title = paste0("GO BP: ", gsub("_", " ", grp)))
      print(p_bp)
      ggsave(file.path(outdir, paste0("GO_BP_", grp, ".png")),
             p_bp, width = 11, height = 8)
    }, error = function(e) {
      cat("    [!] GO BP 可视化失败:", conditionMessage(e), "\n")
    })
    write.csv(go_bp@result,
              file.path(outdir, paste0("GO_BP_", grp, ".csv")),
              row.names = FALSE)
  }

  # GO CC
  go_cc <- do_go(entrez, "CC")
  if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
    cat("    GO CC:", nrow(go_cc@result), "terms\n")
    tryCatch({
      p_cc <- dotplot(go_cc, showCategory = min(10, nrow(go_cc@result)),
                      title = paste0("GO CC: ", gsub("_", " ", grp)))
      print(p_cc)
      ggsave(file.path(outdir, paste0("GO_CC_", grp, ".png")),
             p_cc, width = 11, height = 8)
    }, error = function(e) {
      cat("    [!] GO CC 可视化失败:", conditionMessage(e), "\n")
    })
    write.csv(go_cc@result,
              file.path(outdir, paste0("GO_CC_", grp, ".csv")),
              row.names = FALSE)
  }

  # GO MF
  go_mf <- do_go(entrez, "MF")
  if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
    cat("    GO MF:", nrow(go_mf@result), "terms\n")
    tryCatch({
      p_mf <- dotplot(go_mf, showCategory = min(10, nrow(go_mf@result)),
                      title = paste0("GO MF: ", gsub("_", " ", grp)))
      print(p_mf)
      ggsave(file.path(outdir, paste0("GO_MF_", grp, ".png")),
             p_mf, width = 11, height = 8)
    }, error = function(e) {
      cat("    [!] GO MF 可视化失败:", conditionMessage(e), "\n")
    })
    write.csv(go_mf@result,
              file.path(outdir, paste0("GO_MF_", grp, ".csv")),
              row.names = FALSE)
  }

  # KEGG
  kk <- do_kegg(entrez)
  if (!is.null(kk) && nrow(kk@result) > 0) {
    cat("    KEGG:", nrow(kk@result), "pathways\n")
    tryCatch({
      p_kk <- dotplot(kk, showCategory = min(15, nrow(kk@result)),
                      title = paste0("KEGG: ", gsub("_", " ", grp)))
      print(p_kk)
      ggsave(file.path(outdir, paste0("KEGG_", grp, ".png")),
             p_kk, width = 11, height = 8)
    }, error = function(e) {
      cat("    [!] KEGG 可视化失败:", conditionMessage(e), "\n")
    })
    write.csv(kk@result,
              file.path(outdir, paste0("KEGG_", grp, ".csv")),
              row.names = FALSE)
  }

  # cnetplot
  if (!is.null(go_bp) && nrow(go_bp@result) >= 3) {
    tryCatch({
      # 构建 geneList
      if (grp == "Progressive_Up") {
        gl <- prog_up_table$logFC_NC
        names(gl) <- prog_up_table$gene
      } else {
        gl <- prog_down_table$logFC_NC
        names(gl) <- prog_down_table$gene
      }
      gl <- sort(gl, decreasing = TRUE)
      # 转成 ENTREZID named vector
      s2e_map <- bitr(names(gl), fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
      gl_entrez <- gl[s2e_map$SYMBOL]
      names(gl_entrez) <- s2e_map$ENTREZID

      p_cnet <- cnetplot(go_bp, showCategory = min(5, nrow(go_bp@result)),
                         color.params = list(foldChange = gl_entrez, edge = TRUE),
                         cex.params = list(category_label = 0.8, gene_label = 0.5),
                         circular = FALSE,
                         node_label = "all")
      print(p_cnet)
      ggsave(file.path(outdir, paste0("cnetplot_", grp, ".png")),
             p_cnet, width = 14, height = 10)
      cat("    cnetplot 已保存\n")
    }, error = function(e) {
      cat("    [!] cnetplot 失败:", conditionMessage(e), "\n")
    })
  }

  cat("\n")
}

# ---- GO 比较：Progressive_Up vs Progressive_Down ----
cat("  绘制 GO BP 比较图 ...\n")
if (length(prog_up) >= 5 && length(prog_down) >= 5) {
  comp_list <- list(
    `Progressive Up`   = convert_id(prog_up),
    `Progressive Down` = convert_id(prog_down)
  )
  if (all(sapply(comp_list, length) >= 5)) {
    comp_go <- tryCatch(
      compareCluster(geneClusters = comp_list,
                     fun = "enrichGO", OrgDb = org.Hs.eg.db,
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                     readable = TRUE),
      error = function(e) { cat("    compareCluster失败:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(comp_go) && nrow(comp_go@compareClusterResult) > 0) {
      tryCatch({
        p_comp <- dotplot(comp_go, showCategory = min(12, nrow(comp_go@compareClusterResult)),
                          title = "GO BP: Progressive Up vs Down")
        print(p_comp)
        ggsave(file.path(outdir, "GO_BP_comparison.png"),
               p_comp, width = 12, height = 9)
        cat("  GO_BP_comparison.png 已保存\n")
      }, error = function(e) {
        cat("    [!] compareCluster dotplot 失败:", conditionMessage(e), "\n")
      })
    }
  }
}

# ============================================================
# 6. 已知 CRC 驱动基因的跨阶段表达
# ============================================================
cat("\n>>> 6. 已知 CRC 驱动基因验证 ...\n")

known_drivers <- c(
  "APC", "TP53", "KRAS", "SMAD4", "PIK3CA", "BRAF",
  "MYC", "CTNNB1", "TCF7L2", "AXIN2", "PTEN", "EGFR",
  "VEGFA", "MMP7", "CD44", "CCND1", "MYB", "SOX9",
  "BMP4", "TGFB1", "WNT5A", "FZD7", "LGR5"
)

drivers_in_data <- known_drivers[known_drivers %in% rownames(log2_fpkm)]
drivers_prog_up <- intersect(drivers_in_data, prog_up)
drivers_prog_down <- intersect(drivers_in_data, prog_down)

cat("  数据中检测到:", length(drivers_in_data), "/", length(known_drivers), "\n")
cat("  其中全程上调:", length(drivers_prog_up), ":", paste(drivers_prog_up, collapse = ", "), "\n")
cat("  其中全程下调:", length(drivers_prog_down), ":", paste(drivers_prog_down, collapse = ", "), "\n")

if (length(drivers_in_data) > 0) {
  driver_means <- data.frame()
  for (g in drivers_in_data) {
    vals <- log2_fpkm[g, ]
    driver_means <- rbind(driver_means, data.frame(
      gene    = g,
      Normal  = mean(as.numeric(vals[sample_info$group == "Normal"]),  na.rm = TRUE),
      Adenoma = mean(as.numeric(vals[sample_info$group == "Adenoma"]), na.rm = TRUE),
      CRC     = mean(as.numeric(vals[sample_info$group == "CRC"]),     na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }

  driver_long <- driver_means %>%
    pivot_longer(cols = c(Normal, Adenoma, CRC),
                 names_to = "Stage", values_to = "Expression") %>%
    mutate(Stage = factor(Stage, levels = c("Normal", "Adenoma", "CRC")))

  # 标注是否为全程变化基因
  driver_long$progressive <- ifelse(
    driver_long$gene %in% c(drivers_prog_up, drivers_prog_down),
    "Progressive", "Other"
  )

  p_driver <- ggplot(driver_long, aes(x = Stage, y = Expression, group = gene)) +
    geom_line(aes(color = gene, linewidth = progressive), alpha = 0.8) +
    geom_point(aes(color = gene, size = progressive)) +
    scale_linewidth_manual(values = c("Progressive" = 1.2, "Other" = 0.5), guide = "none") +
    scale_size_manual(values = c("Progressive" = 3, "Other" = 1.5), guide = "none") +
    geom_label_repel(
      data = driver_long %>% filter(Stage == "CRC"),
      aes(label = gene, color = gene),
      size = 3.5, direction = "y", nudge_x = 0.3,
      segment.size = 0.3, max.overlaps = 30, show.legend = FALSE,
      box.padding = 0.3
    ) +
    labs(title = "Known CRC Driver Genes Across Stages",
         subtitle = paste0(length(drivers_in_data), " genes detected; ",
                          length(drivers_prog_up) + length(drivers_prog_down),
                          " are progressive"),
         x = "", y = "Mean log2(FPKM + 1)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

  print(p_driver)
  ggsave(file.path(outdir, "known_drivers_lines.png"),
         p_driver, width = 12, height = 8)
  cat("  known_drivers_lines.png 已保存\n")
}

# ============================================================
# 7. 保存全部结果
# ============================================================
cat("\n>>> 7. 保存结果 ...\n")

save(prog_up, prog_down, prog_up_table, prog_down_table, prog_all,
     na_up, na_down, ac_up, ac_down,
     file = file.path(outdir, "progressive_analysis_results.Rdata"))
cat("  progressive_analysis_results.Rdata 已保存\n")

# ============================================================
# 输出文件清单 + 关键发现
# ============================================================
cat("\n============================================================\n")
cat("  全程上调/下调基因分析 完成！\n")
cat("  输出目录: ", outdir, "\n")
cat("============================================================\n\n")

cat("  === 核心基因表 ===\n")
cat("  progressive_genes.csv\n\n")

cat("  === 代表性拆线图 ===\n")
cat("  progressive_up_lines.png          # 全程上调 top 15\n")
cat("  progressive_down_lines.png        # 全程下调 top 15\n\n")

cat("  === 热图 ===\n")
cat("  progressive_heatmap.png           # 全部全程基因\n\n")

cat("  === 重叠分析 ===\n")
cat("  venn_progressive_up.png           # 上调 Venn\n")
cat("  venn_progressive_down.png         # 下调 Venn\n\n")

cat("  === 富集分析 ===\n")
cat("  GO_BP_Progressive_Up.png / .csv\n")
cat("  GO_CC_Progressive_Up.png / .csv\n")
cat("  GO_MF_Progressive_Up.png / .csv\n")
cat("  KEGG_Progressive_Up.png / .csv\n")
cat("  GO_BP_Progressive_Down.png / .csv\n")
cat("  GO_CC_Progressive_Down.png / .csv\n")
cat("  GO_MF_Progressive_Down.png / .csv\n")
cat("  KEGG_Progressive_Down.png / .csv\n")
cat("  GO_BP_comparison.png\n")
cat("  cnetplot_Progressive_Up.png\n")
cat("  cnetplot_Progressive_Down.png\n\n")

cat("  === 已知驱动基因 ===\n")
cat("  known_drivers_lines.png\n\n")

cat("  === R 数据 ===\n")
cat("  progressive_analysis_results.Rdata\n\n")

# 关键发现
cat("============================================================\n")
cat("  关键发现\n")
cat("============================================================\n")
cat("  全程上调基因:   ", length(prog_up),
    " (在 N→A 和 A→CRC 两个阶段均显著上调)\n")
cat("  全程下调基因:   ", length(prog_down),
    " (在 N→A 和 A→CRC 两个阶段均显著下调)\n")
cat("  合计:           ", length(prog_up) + length(prog_down), "\n\n")
cat("  与 CRC_vs_Normal 的重叠:\n")
cat("    全程上调 ∩ N→CRC 上调: ", length(intersect(prog_up, nc_up)), "\n")
cat("    全程下调 ∩ N→CRC 下调: ", length(intersect(prog_down, nc_down)), "\n")
cat("\n")
