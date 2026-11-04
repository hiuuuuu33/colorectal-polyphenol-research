# Plot Figures 3-4 and supplementary enrichment panels from archived GO/KEGG result tables.

run_enrichment_figures <- function(root = ".", verify = TRUE) {
  required <- c("ggplot2", "patchwork")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
  root <- normalizePath(root, mustWork = TRUE)
  input <- file.path(root, "data", "enrichment", "frozen")
  out <- file.path(root, "results", "enrichment")
  figs <- file.path(root, "figures", "reproduced", "transcriptomics")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(figs, recursive = TRUE, showWarnings = FALSE)
  contrasts <- c("Adenoma_vs_Normal", "CRC_vs_Adenoma", "CRC_vs_Normal")
  categories <- c("GO_BP", "GO_CC", "GO_MF", "KEGG_up", "KEGG_down")
  all_results <- list()
  summary <- data.frame()
  manifests <- data.frame()
  for (contrast in contrasts) {
    plots <- list()
    for (category in categories) {
      filename <- paste0("GSE228012_", category, "_", contrast, ".csv")
      dat <- read.csv(file.path(input, filename), stringsAsFactors = FALSE)
      stopifnot(all(c("ID", "Description", "GeneRatio", "p.adjust", "Count") %in% names(dat)),
                !anyNA(dat$p.adjust))
      significant <- dat[dat$p.adjust < .05, , drop = FALSE]
      write.csv(significant, file.path(out, paste0(category, "_", contrast, "_significant.csv")),
                row.names = FALSE)
      key <- paste(contrast, category, sep = "__")
      all_results[[key]] <- dat
      summary <- rbind(summary, data.frame(contrast = contrast, category = category,
                                           tested_terms = nrow(dat), significant_terms = nrow(significant),
                                           source_mode = "frozen_original_results"))
      manifests <- rbind(manifests, data.frame(input_file = paste0("data/enrichment/frozen/", filename),
                                               records = nrow(dat), mode = "original exported full @result table"))
      shown <- head(significant[order(significant$p.adjust, significant$pvalue), ], 15)
      ratios <- strsplit(shown$GeneRatio, "/", fixed = TRUE)
      shown$GeneRatio_numeric <- vapply(ratios, function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
      shown <- shown[order(shown$GeneRatio_numeric), ]
      shown$label <- vapply(shown$Description, function(x) paste(strwrap(x, width = 42), collapse = "\n"), character(1))
      shown$label <- factor(shown$label, levels = unique(shown$label))
      p <- ggplot2::ggplot(shown, ggplot2::aes(GeneRatio_numeric, label)) +
        ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
        ggplot2::scale_color_gradient(low = "#CB2525", high = "#2866B0", name = "Adjusted P") +
        ggplot2::scale_size_continuous(range = c(3, 9), name = "Gene count") +
        ggplot2::labs(title = gsub("_", " ", category), subtitle = gsub("_", " ", contrast),
                       x = "Gene ratio", y = NULL) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                       axis.text.y = ggplot2::element_text(size = 9),
                       panel.grid.minor = ggplot2::element_blank())
      plots[[category]] <- p
      fig_number <- c(Adenoma_vs_Normal = "03", CRC_vs_Adenoma = "04", CRC_vs_Normal = "S01")[[contrast]]
      panel <- LETTERS[match(category, categories)]
      ggplot2::ggsave(file.path(figs, paste0("figure_", fig_number, panel, "_", category, ".png")),
                      p, width = 10, height = max(4.5, 1.7 + nrow(shown) * .36), dpi = 160,
                      bg = "white")
    }
    combined <- patchwork::wrap_plots(plots, ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
    fig_number <- c(Adenoma_vs_Normal = "03", CRC_vs_Adenoma = "04", CRC_vs_Normal = "S01")[[contrast]]
    ggplot2::ggsave(file.path(figs, paste0("figure_", fig_number, "_enrichment_", contrast, ".png")),
                    combined, width = 20, height = 21, dpi = 140, bg = "white", limitsize = FALSE)
  }
  if (verify) {
    expected_tested <- c(4540L, 381L, 766L, 211L, 283L,
                          3394L, 322L, 563L, 180L, 168L,
                          5630L, 604L, 1050L, 290L, 305L)
    expected_significant <- c(239L, 21L, 68L, 10L, 7L,
                               93L, 19L, 13L, 14L, 1L,
                               711L, 76L, 69L, 19L, 29L)
    stopifnot(identical(summary$tested_terms, expected_tested),
              identical(summary$significant_terms, expected_significant))
    write.csv(data.frame(check = paste(summary$contrast, summary$category, sep = "__"),
                          passed = TRUE, expected_tested = expected_tested,
                          expected_significant = expected_significant),
                file.path(out, "verification.csv"), row.names = FALSE)
  }
  write.csv(summary, file.path(out, "enrichment_summary.csv"), row.names = FALSE)
  write.csv(manifests, file.path(out, "input_manifest.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  message("Enrichment figures complete: 15 frozen tables read; submitted significant-term counts verified.")
  invisible(list(summary = summary, results = all_results))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_enrichment_figures(if (length(args)) args[1] else ".")
}
