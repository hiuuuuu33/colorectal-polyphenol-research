# Read the archived docking workbook and search grids; save Table 5 and score-comparison plots.

run_docking <- function(root = ".") {
  root <- normalizePath(root, mustWork = TRUE)
  for (package in c("readxl", "ggplot2")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Required R package is missing: ", package)
    }
  }
  record_dir <- file.path(root, "data", "docking", "records")
  result_dir <- file.path(root, "results", "docking")
  figure_dir <- file.path(root, "figures", "reproduced", "docking")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  workbook <- file.path(record_dir, "docking_results.xlsx")
  records <- as.data.frame(readxl::read_excel(
    workbook, sheet = "Sheet1", range = "A2:C12", col_names = FALSE,
    col_types = c("text", "numeric", "numeric"), .name_repair = "minimal"
  ))
  names(records) <- c("archived_pose_label", "ADH1C_1HT0", "CXCL8_Fab_6WZM")
  if (nrow(records) != 11L || anyNA(records)) {
    stop("The workbook must contain 11 complete records in Sheet1!A2:C12.")
  }
  records$pubchem_cid <- sub(".*PubChemCID_([0-9]+)_.*", "\\1",
                            records$archived_pose_label)
  if (any(!grepl("^[0-9]+$", records$pubchem_cid)) ||
      anyDuplicated(records$pubchem_cid)) {
    stop("Workbook record labels contain missing or duplicate PubChem CIDs.")
  }
  labels <- c(
    `5280863` = "Kaempferol", `9958` = "4-Methylcatechol",
    `5280343` = "Quercetin", `16129869` = "Punicalagin",
    `5281426` = "Umbelliferone", `5281416` = "Esculetin",
    `19` = "2,3-Dihydroxybenzoic acid", `445154` = "Trans-resveratrol",
    `8530` = "1,4-Naphthoquinone", `4788` = "Phloretin",
    `44584733` = "Punicalagin (supplementary structure)"
  )
  if (!setequal(records$pubchem_cid, names(labels))) {
    stop("Workbook CIDs differ from the submitted Table 5 compound records.")
  }
  records <- records[match(names(labels), records$pubchem_cid), ]
  rownames(records) <- NULL
  records$compound <- unname(labels[records$pubchem_cid])
  records$record_type <- ifelse(records$pubchem_cid == "44584733",
                                "supplementary_structure", "primary")
  records$score_unit <- "kcal/mol"
  records$analysis_status <- "archived_best_mode_score"
  records <- records[, c("compound", "pubchem_cid", "record_type",
                         "ADH1C_1HT0", "CXCL8_Fab_6WZM", "score_unit",
                         "analysis_status", "archived_pose_label")]
  table_path <- file.path(result_dir, "table5_docking_scores.csv")
  utils::write.csv(records, table_path, row.names = FALSE, na = "")

  read_grid <- function(filename, receptor) {
    parameters <- readLines(filename, warn = FALSE)
    parameters <- trimws(parameters[nzchar(trimws(parameters))])
    fields <- strsplit(parameters, "=", fixed = TRUE)
    if (any(lengths(fields) != 2L)) stop("Malformed grid file: ", filename)
    grid <- data.frame(
      receptor = receptor,
      parameter = vapply(fields, function(x) trimws(x[1]), character(1)),
      value = vapply(fields, function(x) as.numeric(trimws(x[2])), numeric(1)),
      stringsAsFactors = FALSE
    )
    expected <- c("center_x", "center_y", "center_z", "size_x", "size_y",
                  "size_z", "exhaustiveness", "num_modes", "energy_range")
    if (!setequal(grid$parameter, expected) || anyNA(grid$value)) {
      stop("Grid parameters differ from the archived nine-parameter format.")
    }
    grid
  }
  grids <- rbind(
    read_grid(file.path(root, "data", "docking", "grids", "ADH1C_1HT0.txt"),
              "ADH1C / 1HT0"),
    read_grid(file.path(root, "data", "docking", "grids", "CXCL8_6WZM.txt"),
              "CXCL8-Fab assembly / 6WZM")
  )
  grid_path <- file.path(result_dir, "archived_grid_parameters.csv")
  utils::write.csv(grids, grid_path, row.names = FALSE)

  long <- rbind(
    data.frame(records, receptor = "ADH1C / 1HT0", score = records$ADH1C_1HT0),
    data.frame(records, receptor = "CXCL8-Fab assembly / 6WZM",
               score = records$CXCL8_Fab_6WZM)
  )
  long$compound <- factor(long$compound, levels = rev(unname(labels)))
  long$receptor <- factor(long$receptor,
                          levels = c("ADH1C / 1HT0", "CXCL8-Fab assembly / 6WZM"))
  score_plot <- ggplot2::ggplot(long, ggplot2::aes(x = score, y = compound)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = score, yend = compound),
                          colour = "grey80", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = receptor, shape = record_type), size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", score)),
                      nudge_y = 0.25, size = 3, show.legend = FALSE) +
    ggplot2::facet_grid(. ~ receptor) +
    ggplot2::scale_colour_manual(values = c("#087F8C", "#B74065"), guide = "none") +
    ggplot2::scale_shape_manual(
      values = c(primary = 16, supplementary_structure = 17),
      breaks = c("primary", "supplementary_structure"),
      labels = c("Primary compound record", "Supplementary punicalagin structure"),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(limits = c(-9, 2.5), breaks = seq(-8, 2, 2)) +
    ggplot2::labs(
      title = "Archived molecular-docking scores",
      subtitle = "Ten prioritized compounds and one supplementary structure record",
      x = "Best-mode AutoDock Vina score (kcal/mol)", y = NULL,
      caption = "Source: docking_results.xlsx. Scores are compared within each receptor system."
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
      legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )
  figure_paths <- file.path(figure_dir, c("docking_score_summary.png", "docking_score_summary.pdf"))
  ggplot2::ggsave(figure_paths[1], score_plot, width = 11.5, height = 7.5,
                  units = "in", dpi = 200, bg = "white")
  ggplot2::ggsave(figure_paths[2], score_plot, width = 11.5, height = 7.5,
                  units = "in", bg = "white")

  status <- data.frame(
    item = c("Table 5 scores", "Grid parameters", "Submitted Figure 9",
             "Vina search", "Pose interaction diagrams"),
    available_inputs = c("Original XLSX workbook", "Original parameter TXT files",
                         "Final submitted PDF", "Raw PDB structures and grid files",
                         "Submitted PDF and archived visualization exports"),
    reproduction = c("Regenerated from workbook", "Parsed into a machine-readable table",
                     "Preserved as a submitted figure", "Requires prepared PDBQT inputs",
                     "Requires final pose coordinates and visualization settings"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(status, file.path(result_dir, "reproduction_status.csv"), row.names = FALSE)
  capture.output(utils::sessionInfo(), file = file.path(result_dir, "sessionInfo.txt"))
  message("Docking records rebuilt: 10 primary records + 1 supplementary structure, 22 scores.")
  invisible(list(scores = records, grids = grids, status = status,
                 table_paths = c(table_path, grid_path), figure_paths = figure_paths,
                 plot = score_plot))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_docking(if (length(args)) args[1] else ".")
}
