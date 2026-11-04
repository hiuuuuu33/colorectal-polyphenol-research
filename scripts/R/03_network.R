# Redraw manuscript Figures 5, 6 and reconstructed Figure 8 from frozen inputs.
# Run network_analysis.py once before this script if results are absent.

run_network_figures <- function(project_root = ".") {
  packages <- c("ggplot2", "ggrepel")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install required R packages: ", paste(absent, collapse = ", "))
  library(ggplot2)
  root <- normalizePath(project_root, mustWork = TRUE)
  inputs <- file.path(root, "data", "network")
  results <- file.path(root, "results", "network")
  figures <- file.path(root, "figures", "reproduced", "network")
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  load_csv <- function(name) read.csv(file.path(results, name), stringsAsFactors = FALSE)
  save_figure <- function(plot, name, width, height) {
    ggsave(file.path(figures, paste0(name, ".png")), plot,
           width = width, height = height, dpi = 220, bg = "white")
    ggsave(file.path(figures, paste0(name, ".pdf")), plot,
           width = width, height = height, device = grDevices::pdf)
  }
  theme_set(theme_bw(base_size = 11) +
              theme(panel.grid.minor = element_blank(),
                    plot.title = element_text(face = "bold", size = 13)))

  nodes <- load_csv("ppi_centrality_recalculated.csv")
  positions <- load_csv("ppi_layout.csv")
  nodes <- merge(nodes, positions, by = "name", sort = FALSE)
  edges <- read.delim(file.path(inputs, "string_interactions_short.tsv"),
                      check.names = FALSE)
  a <- match(edges[[1]], nodes$name)
  b <- match(edges[[2]], nodes$name)
  segments <- data.frame(x = nodes$x[a], y = nodes$y[a],
                         xend = nodes$x[b], yend = nodes$y[b])
  network <- ggplot() +
    geom_segment(data = segments, aes(x, y, xend = xend, yend = yend),
                 color = "#AAB4B8", linewidth = 0.5, alpha = 0.7) +
    geom_point(data = nodes, aes(x, y, size = degree, fill = degree),
               shape = 21, color = "white", stroke = 0.6) +
    ggrepel::geom_text_repel(data = nodes, aes(x, y, label = name),
                             seed = 2026, size = 3.3, max.overlaps = Inf,
                             min.segment.length = 0, box.padding = 0.3) +
    scale_fill_gradient(low = "#AFD5D0", high = "#157D80") +
    scale_size_continuous(range = c(3, 9), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    labs(title = "Progression-associated protein network", size = "Degree", fill = "Degree") +
    coord_cartesian(clip = "off") + theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13), plot.margin = margin(18, 24, 18, 24))
  save_figure(network, "figure_05A_ppi_network", 11, 6.8)

  top <- nodes[order(-nodes$degree, -nodes$betweenness, nodes$name), ][1:10, ]
  metrics <- rbind(data.frame(name = top$name, metric = "Degree", value = top$degree),
                   data.frame(name = top$name, metric = "Betweenness", value = top$betweenness),
                   data.frame(name = top$name, metric = "Closeness", value = top$closeness))
  metrics$name <- factor(metrics$name, levels = rev(top$name))
  centrality <- ggplot(metrics, aes(value, name, fill = metric)) +
    geom_col(width = 0.7) + facet_wrap(~metric, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = c("Betweenness" = "#8A6792", "Closeness" = "#B49342", "Degree" = "#157D80")) +
    labs(title = "Topology of the ten highest-degree proteins", x = NULL, y = NULL) +
    theme(legend.position = "none", panel.grid.major.y = element_blank())
  save_figure(centrality, "figure_05B_ppi_centrality", 11, 5.2)

  all_records <- load_csv("proximity_all_66_records.csv")
  all_records <- all_records[order(all_records$z_closest), ]
  all_records$record_rank <- seq_len(nrow(all_records))
  distribution <- ggplot(all_records, aes(record_rank, z_closest)) +
    geom_hline(yintercept = -1.96, color = "#B34C50", linetype = "dashed") +
    geom_hline(yintercept = 1.96, color = "#B34C50", linetype = "dashed") +
    geom_point(aes(fill = z_closest), shape = 21, size = 3.2, color = "white", stroke = 0.4) +
    scale_fill_gradient2(low = "#157D80", mid = "#E5E8E6", high = "#B34C50", midpoint = 0) +
    labs(title = "Network proximity across 66 recorded entries", x = "Entry ordered by closest-distance Z score",
         y = expression(Z[closest])) + theme(legend.position = "none")
  save_figure(distribution, "figure_06A_proximity_distribution", 10, 4.6)

  top10 <- load_csv("top10_candidates.csv")
  top10$display_name <- factor(top10$display_name, levels = rev(top10$display_name))
  lollipop <- ggplot(top10, aes(-z_closest, display_name)) +
    geom_segment(aes(xend = 0, yend = display_name), linewidth = 0.9, color = "#B7C4C3") +
    geom_point(aes(size = n_mapped_chemical, fill = -z_closest), shape = 21, color = "white", stroke = 0.8) +
    geom_text(aes(label = sprintf("%.2f", z_closest)), nudge_x = 0.2, hjust = 0, size = 3.5) +
    scale_size_continuous(range = c(3.5, 10), name = "Mapped targets",
                          guide = guide_legend(override.aes = list(fill = "#157D80", color = "#157D80"))) +
    scale_fill_gradient(low = "#B9DAD5", high = "#157D80", name = expression(-Z[closest])) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
    labs(title = "Ten prioritized compound identities", x = expression(-Z[closest]), y = NULL) +
    theme(panel.grid.major.y = element_blank())
  save_figure(lollipop, "figure_06B_top10_candidates", 11, 6.5)

  string_layout <- file.path(results, "compound_string_layouts.csv")
  string_summary <- file.path(results, "compound_string_summary.csv")
  if (file.exists(string_layout) && file.exists(string_summary)) {
    layouts <- read.csv(string_layout, stringsAsFactors = FALSE)
    compound_info <- read.csv(string_summary, stringsAsFactors = FALSE)
    compound_edges <- list()
    for (i in seq_len(nrow(compound_info))) {
      item <- compound_info[i, ]
      dat <- read.delim(file.path(root, item$edges_file), check.names = FALSE)
      pos <- layouts[layouts$chemical == item$chemical, ]
      dat <- dat[dat$stringId_A %in% pos$string_id & dat$stringId_B %in% pos$string_id, ]
      if (nrow(dat)) {
        ai <- match(dat$stringId_A, pos$string_id)
        bi <- match(dat$stringId_B, pos$string_id)
        compound_edges[[i]] <- data.frame(chemical = item$chemical,
                                          x = pos$x[ai], y = pos$y[ai],
                                          xend = pos$x[bi], yend = pos$y[bi])
      }
    }
    segment_data <- do.call(rbind, compound_edges)
    labels <- setNames(paste0(LETTERS[seq_len(nrow(compound_info))], ". ",
                              compound_info$display_name), compound_info$chemical)
    layouts$chemical <- factor(layouts$chemical, levels = compound_info$chemical)
    segment_data$chemical <- factor(segment_data$chemical, levels = compound_info$chemical)
    compound_plot <- ggplot() +
      geom_segment(data = segment_data, aes(x, y, xend = xend, yend = yend),
                   color = "#BAC5C8", linewidth = 0.22, alpha = 0.55) +
      geom_point(data = layouts, aes(x, y, size = degree), color = "#157D80", alpha = 0.9) +
      geom_text(data = layouts, aes(x, y, label = symbol), size = 1.7, vjust = -0.55,
                check_overlap = TRUE) +
      facet_wrap(~chemical, ncol = 2, labeller = labeller(chemical = labels)) +
      scale_size_continuous(range = c(1, 3)) +
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.18)) +
      labs(title = "Recorded compound targets: STRING 12.0 functional associations") +
      theme_void(base_size = 11) +
      theme(legend.position = "none", strip.text = element_text(hjust = 0, face = "bold"),
             panel.spacing = grid::unit(1.2, "lines"),
             plot.title = element_text(size = 13, face = "bold", margin = margin(b = 15)))
    save_figure(compound_plot, "figure_08_compound_target_networks_reconstructed", 11, 16)
  }
  message("Network figures saved to ", figures)
  invisible(list(nodes = nodes, candidates = top10))
}

if (sys.nframe() == 0L) run_network_figures()
