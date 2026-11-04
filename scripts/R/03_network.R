# Display native STRING networks and plot protein centrality and compound proximity for Figures 5, 6, and 8.

run_network_figures <- function(project_root = ".") {
  packages <- c("ggplot2")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install required R packages: ", paste(absent, collapse = ", "))
  library(ggplot2)
  root <- normalizePath(project_root, mustWork = TRUE)
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
  python <- Sys.getenv("PYTHON", unset = Sys.which("python3"))
  if (!nzchar(python)) stop("Python 3 is required; set PYTHON to its executable.")
  native_log <- system2(python, c(file.path(root, "scripts", "build_string_visuals.py"), "--offline"),
                         stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(native_log, "status"))) stop(paste(native_log, collapse = "\n"))

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

  message("Network figures saved to ", figures)
  invisible(list(nodes = nodes, candidates = top10))
}

if (sys.nframe() == 0L) run_network_figures()
