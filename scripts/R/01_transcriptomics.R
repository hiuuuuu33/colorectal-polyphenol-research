# Reproduce the submitted analysis from the supplied processed FPKM matrix.
# Run from the repository root: Rscript scripts/R/01_transcriptomics.R

run_transcriptomics <- function(root = ".", verify = TRUE) {
  required <- c("limma", "ggplot2", "ggrepel", "pheatmap", "patchwork")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
  root <- normalizePath(root, mustWork = TRUE)
  input_dir <- file.path(root, "data", "transcriptomics")
  out <- file.path(root, "results", "transcriptomics")
  figs <- file.path(root, "figures", "reproduced", "transcriptomics")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(figs, recursive = TRUE, showWarnings = FALSE)
  set.seed(20261104)

  raw <- read.table(file.path(input_dir, "GSE228012_all.gene.FPKM.mRNA_processed.txt"),
                    header = TRUE, sep = "\t", check.names = TRUE)
  symbols <- as.character(raw[[2]])
  fpkm <- as.matrix(raw[, -(1:2)])
  storage.mode(fpkm) <- "double"
  ranked <- order(rowMeans(fpkm), decreasing = TRUE)
  retained <- ranked[!duplicated(symbols[ranked])]
  fpkm <- fpkm[retained, , drop = FALSE]
  rownames(fpkm) <- symbols[retained]
  gene_ids <- data.frame(gene = symbols[retained], GeneID = raw[[1]][retained])
  colnames(fpkm) <- trimws(sub("\\.FPKM$", "", colnames(fpkm)))
  group <- ifelse(grepl("T$", colnames(fpkm)), "CRC",
                  ifelse(grepl("N$", colnames(fpkm)), "Normal",
                         ifelse(grepl("pre$", colnames(fpkm)), "Adenoma", NA_character_)))
  if (anyNA(group)) stop("Unrecognized sample identifier.")
  group <- factor(group, levels = c("Normal", "Adenoma", "CRC"))
  sample_info <- data.frame(sample = colnames(fpkm), group = group,
                            row.names = colnames(fpkm))
  keep <- rowSums(fpkm > 1) >= round(0.2 * ncol(fpkm))
  log2_fpkm <- log2(fpkm[keep, , drop = FALSE] + 1)
  write.csv(sample_info, file.path(out, "sample_metadata.csv"), row.names = FALSE)
  write.csv(gene_ids[keep, ], file.path(out, "tested_gene_identifiers.csv"), row.names = FALSE)

  # This is the original unblocked limma model used for the submitted figures.
  design <- stats::model.matrix(~0 + group)
  colnames(design) <- levels(group)
  contrasts <- limma::makeContrasts(
    Adenoma_vs_Normal = Adenoma - Normal,
    CRC_vs_Adenoma = CRC - Adenoma,
    CRC_vs_Normal = CRC - Normal, levels = design)
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(log2_fpkm, design), contrasts))
  de <- lapply(colnames(contrasts), function(nm) {
    tab <- limma::topTable(fit, coef = nm, number = Inf, adjust.method = "BH")
    tab$gene <- rownames(tab)
    tab <- tab[, c("gene", setdiff(names(tab), "gene"))]
    tab$direction <- ifelse(tab$adj.P.Val < .05 & tab$logFC > 1, "Up",
                             ifelse(tab$adj.P.Val < .05 & tab$logFC < -1, "Down", "Stable"))
    write.csv(tab, file.path(out, paste0("DE_", nm, "_all.csv")), row.names = FALSE)
    write.csv(tab[tab$direction != "Stable", ],
              file.path(out, paste0("DE_", nm, "_significant.csv")), row.names = FALSE)
    tab
  })
  names(de) <- colnames(contrasts)
  summary <- do.call(rbind, lapply(names(de), function(nm) {
    tab <- de[[nm]]
    data.frame(contrast = nm, tested = nrow(tab), significant = sum(tab$direction != "Stable"),
               up = sum(tab$direction == "Up"), down = sum(tab$direction == "Down"))
  }))
  write.csv(summary, file.path(out, "differential_expression_summary.csv"), row.names = FALSE)

  na <- de$Adenoma_vs_Normal
  ac <- de$CRC_vs_Adenoma
  nc <- de$CRC_vs_Normal
  up <- intersect(na$gene[na$direction == "Up"], ac$gene[ac$direction == "Up"])
  down <- intersect(na$gene[na$direction == "Down"], ac$gene[ac$direction == "Down"])
  build_progressive <- function(genes, direction) {
    data.frame(gene = genes,
               logFC_NA = na$logFC[match(genes, na$gene)],
               adjP_NA = na$adj.P.Val[match(genes, na$gene)],
               logFC_AC = ac$logFC[match(genes, ac$gene)],
               adjP_AC = ac$adj.P.Val[match(genes, ac$gene)],
               logFC_NC = nc$logFC[match(genes, nc$gene)],
               adjP_NC = nc$adj.P.Val[match(genes, nc$gene)], direction = direction)
  }
  prog <- rbind(build_progressive(up, "Up-Up"), build_progressive(down, "Down-Down"))
  write.csv(prog, file.path(out, "progressive_genes.csv"), row.names = FALSE)
  preprocessing <- data.frame(step = c("Input rows", "Unique gene symbols", "Tested genes", "Samples",
                                       "Progressive up", "Progressive down"),
                              value = c(nrow(raw), nrow(fpkm), nrow(log2_fpkm), ncol(log2_fpkm),
                                        length(up), length(down)))
  write.csv(preprocessing, file.path(out, "preprocessing_summary.csv"), row.names = FALSE)

  checks <- data.frame(check = character(), passed = logical(), max_absolute_error = numeric())
  if (verify) {
    expected <- data.frame(significant = c(725L, 312L, 1544L), up = c(254L, 178L, 855L),
                           down = c(471L, 134L, 689L))
    stopifnot(nrow(raw) == 18116L, nrow(fpkm) == 18114L, nrow(log2_fpkm) == 12298L,
              all(table(group) == 11L), identical(summary[, c("significant", "up", "down")], expected),
              length(up) == 8L, length(down) == 21L)
    refs <- c(Adenoma_vs_Normal = "adenoma_normal_diff.csv", CRC_vs_Adenoma = "crc_adenoma_diff.csv",
               CRC_vs_Normal = "crc_normal_diff.csv")
    for (nm in names(refs)) {
      ref <- read.csv(file.path(input_dir, "reference", refs[[nm]]), row.names = 1)
      observed <- de[[nm]][de[[nm]]$direction != "Stable", ]
      ix <- match(ref$gene, observed$gene)
      if (anyNA(ix) || nrow(ref) != nrow(observed)) stop("DEG identity mismatch: ", nm)
      err <- max(abs(ref$logFC - observed$logFC[ix]), abs(ref$adj.P.Val - observed$adj.P.Val[ix]))
      stopifnot(err < 1e-10)
      checks <- rbind(checks, data.frame(check = paste0(nm, "_saved_DEG_table"), passed = TRUE,
                                         max_absolute_error = err))
    }
    ref_prog <- read.csv(file.path(input_dir, "reference", "progressive_genes.csv"))
    stopifnot(setequal(ref_prog$gene, prog$gene))
    checks <- rbind(checks, data.frame(check = "29_progressive_gene_identities", passed = TRUE,
                                       max_absolute_error = 0))
    write.csv(checks, file.path(out, "verification.csv"), row.names = FALSE)
  }

  colors <- c(Normal = "#00A087", Adenoma = "#E64B35", CRC = "#3C5488")
  plot_theme <- ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.title = ggplot2::element_blank())
  save_plot <- function(plot, filename, width = 9, height = 7) {
    ggplot2::ggsave(file.path(figs, filename), plot, width = width, height = height, dpi = 160,
                    bg = "white", limitsize = FALSE)
  }

  # FactoMineR::PCA scales each gene to unit variance by default.
  pca <- stats::prcomp(t(log2_fpkm), center = TRUE, scale. = TRUE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pc <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group)
  write.csv(pc, file.path(out, "pca_coordinates.csv"), row.names = FALSE)
  write.csv(data.frame(component = seq_along(variance), variance_percent = variance),
            file.path(out, "pca_variance.csv"), row.names = FALSE)
  pca_plot <- ggplot2::ggplot(pc, ggplot2::aes(PC1, PC2, color = group, shape = group)) +
    ggplot2::stat_ellipse(type = "norm", linewidth = .5, show.legend = FALSE) +
    ggplot2::geom_point(size = 3) + ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(title = "PCA of GSE228012",
                   x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2])) + plot_theme
  save_plot(pca_plot, "figure_01A_PCA.png", 8, 6)
  volcano_plots <- lapply(names(de), function(nm) {
    dat <- de[[nm]]
    dat$neg_log10_adjP <- -log10(pmax(dat$adj.P.Val, .Machine$double.xmin))
    upper <- dat[dat$direction == "Up", ]
    lower <- dat[dat$direction == "Down", ]
    labels <- rbind(head(upper[order(-upper$logFC), ], 5), head(lower[order(lower$logFC), ], 5))
    p <- ggplot2::ggplot(dat, ggplot2::aes(logFC, neg_log10_adjP)) +
      ggplot2::geom_point(ggplot2::aes(color = direction), alpha = .4, size = 1.2) +
      ggplot2::scale_color_manual(values = c(Up = "#DC0000", Down = "#3C5488", Stable = "grey70")) +
      ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 4, color = "grey40", linewidth = .5) +
      ggplot2::geom_hline(yintercept = -log10(.05), linetype = 4, color = "grey40", linewidth = .5) +
      ggrepel::geom_label_repel(data = labels, ggplot2::aes(label = gene),
                                seed = 20261104, color = "black", size = 2.8,
                                max.overlaps = Inf, box.padding = .4) +
      ggplot2::labs(title = gsub("_", " ", nm), x = "log2 fold change", y = "-log10 adjusted P") + plot_theme
    panel <- c(Adenoma_vs_Normal = "B", CRC_vs_Adenoma = "C", CRC_vs_Normal = "D")[[nm]]
    save_plot(p, paste0("figure_01", panel, "_volcano_", nm, ".png"))
    p
  })
  names(volcano_plots) <- names(de)
  figure1 <- patchwork::wrap_plots(c(list(pca_plot), volcano_plots), ncol = 2) +
    patchwork::plot_annotation(tag_levels = "A")
  save_plot(figure1, "figure_01_transcriptomic_overview.png", 16, 12)

  stage_order <- order(group)
  heat <- log2_fpkm[prog$gene, stage_order, drop = FALSE]
  heat_z <- t(scale(t(heat)))
  col_anno <- data.frame(Stage = group[stage_order], row.names = colnames(heat))
  row_anno <- data.frame(Direction = prog$direction, row.names = prog$gene)
  pheat <- pheatmap::pheatmap(heat_z, silent = TRUE, cluster_cols = FALSE, cluster_rows = TRUE,
                             show_colnames = FALSE, fontsize_row = 8, annotation_col = col_anno,
                             annotation_row = row_anno,
                             annotation_colors = list(Stage = colors,
                                                       Direction = c("Up-Up" = "#DC0000", "Down-Down" = "#3C5488")),
                             gaps_col = c(11, 22), main = "29 progression-associated genes")
  png(file.path(figs, "figure_02A_progressive_heatmap.png"), width = 1500, height = 1250, res = 160)
  grid::grid.newpage(); grid::grid.draw(pheat$gtable); grDevices::dev.off()
  write.csv(data.frame(gene = rownames(heat_z), heat_z, check.names = FALSE),
            file.path(out, "progressive_heatmap_zscores.csv"), row.names = FALSE)

  means <- data.frame(gene = prog$gene, direction = prog$direction,
                      Normal = rowMeans(log2_fpkm[prog$gene, group == "Normal", drop = FALSE]),
                      Adenoma = rowMeans(log2_fpkm[prog$gene, group == "Adenoma", drop = FALSE]),
                      CRC = rowMeans(log2_fpkm[prog$gene, group == "CRC", drop = FALSE]))
  write.csv(means, file.path(out, "progressive_stage_means.csv"), row.names = FALSE)
  long <- reshape(means, varying = c("Normal", "Adenoma", "CRC"), v.names = "Expression",
                   timevar = "Stage", times = c("Normal", "Adenoma", "CRC"), direction = "long")
  long$Stage <- factor(long$Stage, levels = levels(group))
  line_plot <- function(dir, color, title) {
    # The submitted downregulated profile displays the 15 largest overall changes.
    chosen <- prog[prog$direction == dir, ]
    chosen <- head(chosen[order(-abs(chosen$logFC_NC)), "gene"], 15)
    dat <- long[long$gene %in% chosen, ]
    ggplot2::ggplot(dat, ggplot2::aes(Stage, Expression, group = gene)) +
      ggplot2::geom_line(color = color, linewidth = .7, alpha = .8) +
      ggplot2::geom_point(color = color, size = 1.8) +
      ggrepel::geom_label_repel(data = dat[dat$Stage == "CRC", ], ggplot2::aes(label = gene),
                                seed = 20261104, direction = "y", nudge_x = .35,
                                max.overlaps = Inf, size = 2.6, color = color) +
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(.2, 1.2))) +
      ggplot2::labs(title = title, x = NULL, y = "Mean log2(FPKM + 1)") + plot_theme
  }
  down_plot <- line_plot("Down-Down", "#3C5488", "Progressively downregulated genes (top 15)")
  up_plot <- line_plot("Up-Up", "#DC0000", "Progressively upregulated genes (all 8)")
  save_plot(down_plot, "figure_02B_progressive_down_lines.png", 10, 7)
  save_plot(up_plot, "figure_02C_progressive_up_lines.png", 10, 7)
  heat_plot <- patchwork::wrap_elements(full = pheat$gtable)
  figure2 <- patchwork::wrap_plots(heat_plot, down_plot, up_plot,
                                  design = "AB\nAC", widths = c(1.15, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_plot(figure2, "figure_02_progressive_genes.png", 16, 11)

  saveRDS(list(log2_fpkm = log2_fpkm, sample_info = sample_info, de = de, progressive = prog,
               pca_variance = variance), file.path(out, "transcriptomics_results.rds"))
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  message("Transcriptomics complete: 725 / 312 / 1544 DEGs; 8 up + 21 down progression-associated genes.")
  invisible(list(summary = summary, progressive = prog, preprocessing = preprocessing, verification = checks))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_transcriptomics(if (length(args)) args[1] else ".")
}
