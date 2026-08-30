#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "PANDA plotting command\n\n",
    "Usage:\n",
    "  Rscript cli/panda_plot.R --results RESULTS_DIR\n",
    "  Rscript cli/panda_plot.R --results RESULTS_DIR --figures heatmap\n",
    "  Rscript cli/panda_plot.R --results RESULTS_DIR --figures asm samples\n",
    "  Plot types: distribution, heatmap, lollipop, asm, samples, group\n"
  )
}

if (!length(args) || "--help" %in% args || "-h" %in% args) {
  usage()
  quit(status = 0L)
}

if (!"--results" %in% args) {
  usage()
  stop("--results is required.", call. = FALSE)
}

results_index <- match("--results", args)
if (results_index == length(args)) {
  stop("--results requires a directory.", call. = FALSE)
}
results_dir <- normalizePath(args[[results_index + 1L]], mustWork = TRUE)
json_path <- file.path(results_dir, "PANDA_analysis.json")
if (!file.exists(json_path)) {
  stop("PANDA_analysis.json was not found in: ", results_dir, call. = FALSE)
}

plot_spec <- "all"
if (any(c("--figure", "--plots") %in% args)) {
  stop("Unknown plotting option. Use --figures followed by one or more plot types.", call. = FALSE)
}
if ("--figures" %in% args) {
  figure_index <- match("--figures", args)
  if (figure_index == length(args)) stop("--figures requires at least one plot type.", call. = FALSE)
  following <- args[(figure_index + 1L):length(args)]
  stop_at <- which(grepl("^-", following))
  if (length(stop_at)) following <- following[seq_len(stop_at[[1L]] - 1L)]
  if (!length(following)) stop("--figures requires at least one plot type.", call. = FALSE)
  plot_spec <- tolower(following)
}

top_n_explicit <- "--top-n" %in% args
if (top_n_explicit && !("all" %in% plot_spec || "lollipop" %in% plot_spec)) {
  stop(
    "--top-n is only applicable when --figures all or --figures lollipop is selected. ",
    "Distribution, heatmap, ASM, samples, and group plots always use all retained variants.",
    call. = FALSE
  )
}

top_n <- 30L
if ("--top-n" %in% args) {
  top_index <- match("--top-n", args)
  if (top_index == length(args)) stop("--top-n requires an integer.", call. = FALSE)
  top_n <- as.integer(args[[top_index + 1L]])
}
if (is.na(top_n) || top_n < 1L) stop("--top-n must be positive.", call. = FALSE)

valid_plots <- c("distribution", "heatmap", "lollipop", "asm", "samples", "group")
plots <- if ("all" %in% plot_spec) valid_plots else unique(plot_spec)
unknown <- setdiff(plots, valid_plots)
if (length(unknown)) {
  stop("Unknown plot type(s): ", paste(unknown, collapse = ", "),
       ". Choose: all, distribution, heatmap, lollipop, asm, samples, group.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## Cairo gives substantially cleaner antialiasing for the very thin
## Count-weighted bands in the abundance heatmap PDF.  Fall back to the
## standard device on systems where Cairo PDF is unavailable.
pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"

bundle <- jsonlite::read_json(json_path, simplifyVector = TRUE)
plot_dir <- file.path(results_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

as_data_frame <- function(x) {
  if (is.null(x)) return(data.frame())
  if (is.data.frame(x)) return(x)
  if (!length(x)) return(data.frame())
  bind_rows(lapply(x, function(row) as.data.frame(row, stringsAsFactors = FALSE)))
}

sample_records <- bundle$samples
if (is.null(names(sample_records))) stop("No sample records found in JSON.", call. = FALSE)

for (sample_id in names(sample_records)) {
  rec <- sample_records[[sample_id]]
  long_data <- as_data_frame(rec$long_data)
  if (!nrow(long_data)) next
  if (!"Count" %in% names(long_data)) long_data$Count <- 1L
  long_data$Count[is.na(long_data$Count) | long_data$Count < 1] <- 1L
  
  ## Count is the objective criterion for selecting Top-N display rows.
  ## Mean methylation is used only as a deterministic tie-breaker/order.
  counts <- long_data %>%
    group_by(ReadID) %>%
    summarise(
      Count = first(Count),
      Mean_Meth = mean(Methylation, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Count), desc(Mean_Meth), ReadID) %>%
    mutate(Count_Rank = row_number())
  ## Distribution and heatmap use all retained variants. Top-N is only a
  ## readability control for the row-oriented lollipop plot.
  plot_data <- long_data
  plot_counts <- counts
  lollipop_ids <- counts %>% slice_head(n = top_n) %>% pull(ReadID)
  lollipop_data <- long_data %>% filter(ReadID %in% lollipop_ids)
  ## Reorder selected rows for visual interpretation, without changing
  ## which rows were selected.
  lollipop_counts <- counts %>% filter(ReadID %in% lollipop_ids) %>%
    arrange(Mean_Meth, desc(Count), ReadID)
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
  
  read_stats <- plot_data %>%
    group_by(ReadID) %>%
    summarise(Mean_Meth = mean(Methylation, na.rm = TRUE) * 100,
              Count = first(Count), .groups = "drop") %>%
    ## Invalid/empty methylation patterns should not reach geom_histogram;
    ## filtering here avoids misleading per-sample geom_bar warnings.
    filter(is.finite(Mean_Meth), Mean_Meth >= 0, Mean_Meth <= 100,
           is.finite(Count), Count > 0)
  
  if ("distribution" %in% plots) {
    p <- ggplot(read_stats, aes(x = Mean_Meth, weight = Count)) +
      geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
      scale_x_continuous(breaks = seq(0, 100, 10)) +
      coord_cartesian(xlim = c(0, 100)) +
      theme_minimal(base_size = 13) +
      labs(title = paste0(sample_id, ": methylation distribution"),
           subtitle = "All retained variants; bars are weighted by Count",
           x = "Mean methylation per variant (%)", y = "Weighted read count")
    ggsave(file.path(plot_dir, paste0(safe_name, "_methylation_distribution.pdf")),
           p, width = 8, height = 6, device = pdf_device)
  }
  
  if ("heatmap" %in% plots) {
    ## Preserve the original abundance encoding: each epiallele is a band
    ## whose thickness is proportional to Count, with a fixed page height.
    heat_order <- counts %>% arrange(Mean_Meth, desc(Count), ReadID) %>%
      mutate(ymax = cumsum(Count), ymin = lag(ymax, default = 0)) %>%
      select(ReadID, Count, Mean_Meth, ymin, ymax)
    heat_data <- plot_data %>% select(-any_of("Count")) %>%
      left_join(heat_order, by = "ReadID")
    heat_data$Position <- as.numeric(as.character(heat_data$Position))
    cpg_positions <- sort(unique(heat_data$Position[is.finite(heat_data$Position)]))
    x_width <- if (length(cpg_positions) > 1L) min(5, 0.8 * min(diff(cpg_positions))) else 5
    heat_data$Mean_Meth_Pct <- 100 * heat_order$Mean_Meth[
      match(heat_data$ReadID, heat_order$ReadID)
    ]
    meth_levels <- sprintf("%d-%d%%", seq(0, 90, 10), seq(10, 100, 10))
    heat_data$Meth_Bin <- cut(heat_data$Mean_Meth_Pct,
                              breaks = seq(0, 100, 10), include.lowest = TRUE,
                              labels = meth_levels, right = TRUE)
    binned <- heat_data %>% group_by(Meth_Bin, Position) %>%
      summarise(Weighted_Methylation = weighted.mean(Methylation, w = Count, na.rm = TRUE),
                .groups = "drop")
    binned$Meth_Bin <- factor(binned$Meth_Bin, levels = meth_levels)
    p <- ggplot(binned, aes(x = Position, y = Meth_Bin, fill = Weighted_Methylation * 100)) +
      geom_tile(width = x_width, height = 1) +
      scale_fill_gradient2(low = "lightblue", mid = "white", high = "firebrick",
                           midpoint = 50, limits = c(0, 100), name = "Methylation (%)") +
      scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
      labs(title = paste0(sample_id, ": binned abundance heatmap"),
           subtitle = "All retained reads; display-only aggregation by mean methylation",
           x = "CpG position (bp)", y = "Mean methylation per epiallele")
    heat_height <- 7
    p <- p + theme_minimal(base_size = 18) +
      theme(panel.background = element_rect(fill = "grey92", colour = NA),
            panel.grid.major = element_line(colour = "white", linewidth = 0.35),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 17))
    ggsave(file.path(plot_dir, paste0(safe_name, "_abundance_heatmap.pdf")),
           p, width = 13, height = heat_height,
           device = pdf_device)
    
  }
  
  if ("lollipop" %in% plots) {
    lollipop_data$CpG_Index <- match(lollipop_data$Position, sort(unique(lollipop_data$Position)))
    lollipop_data$ReadID <- factor(lollipop_data$ReadID, levels = lollipop_counts$ReadID)
    lollipop_counts$ReadID <- factor(lollipop_counts$ReadID,
                                     levels = levels(lollipop_data$ReadID))
    plot_title <- paste0(sample_id, ": methylation lollipop (Top ",
                         min(top_n, nrow(lollipop_counts)), ")")
    plot_subtitle <- "Top-N selected by Count; rows displayed by methylation state"
    p1 <- ggplot(lollipop_data, aes(x = CpG_Index, y = ReadID)) +
      geom_line(aes(group = ReadID), color = "grey80") +
      geom_point(aes(fill = factor(Methylation)), shape = 21, size = 3, color = "black") +
      scale_fill_manual(values = c("0" = "white", "1" = "black"), guide = "none") +
      theme_minimal(base_size = 16) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14),
            axis.title = element_text(size = 15)) +
      labs(x = "CpG index", y = "Variant")
    ## Restore the abundance panel used in the GUI.  It is informative for
    ## dereplicated amplicon/NGS data; for Sanger clones (Count == 1) it is
    ## omitted because every bar would have identical height.
    is_sanger <- all(lollipop_counts$Count == 1L)
    p <- p1
    if (!is_sanger) {
      p2 <- ggplot(lollipop_counts, aes(x = Count, y = ReadID)) +
        geom_col(fill = "steelblue") +
        geom_text(aes(label = paste0("#", Count_Rank, " (n=", Count, ")")),
                  hjust = -0.1, size = 4.2) +
        scale_x_continuous(expand = expansion(mult = c(0, 1.05))) +
        theme_void(base_size = 13) +
        theme(plot.margin = margin(l = 10, r = 12)) +
        coord_cartesian(clip = "off") +
        labs(title = "Read Count")
      p <- p1 + p2 + plot_layout(widths = c(3, 1.35))
    }
    p <- p + plot_annotation(
      title = plot_title,
      subtitle = plot_subtitle,
      theme = theme(
        plot.title = element_text(size = 18, face = "plain"),
        plot.subtitle = element_text(size = 14)
      )
    )
    ggsave(file.path(plot_dir, paste0(safe_name, "_lollipop.pdf")),
           p, width = if (is_sanger) 10.5 else 14.2, height = 8.3,
           device = pdf_device)
  }

  if ("asm" %in% plots) {
    detail <- rec$heterogeneity_detail
    mat_df <- if (!is.null(detail$meth_mat)) as_data_frame(detail$meth_mat) else data.frame()
    if (nrow(mat_df) > 1L && "ReadID" %in% names(mat_df)) {
      ids <- as.character(mat_df$ReadID)
      mat_df$ReadID <- NULL
      mat <- as.matrix(mat_df)
      storage.mode(mat) <- "numeric"
      cl <- detail$clusters
      if (!is.null(cl)) {
        cl <- unlist(cl, use.names = TRUE)
        if (!is.null(names(cl))) cl <- as.integer(cl[ids]) else cl <- as.integer(cl)
      }
      if (is.null(cl) || length(cl) != length(ids) || anyNA(cl)) cl <- rep(1L, length(ids))
      row_order <- order(cl, rowMeans(mat, na.rm = TRUE))
      long_asm <- as.data.frame(mat[row_order, , drop = FALSE], check.names = FALSE) %>%
        mutate(ReadID = ids[row_order], Cluster = factor(cl[row_order])) %>%
        tidyr::pivot_longer(cols = -c(ReadID, Cluster), names_to = "Position", values_to = "Methylation")
      long_asm$Position <- factor(long_asm$Position, levels = colnames(mat))
      p_asm <- ggplot(long_asm, aes(x = Position, y = ReadID, fill = factor(Methylation))) +
        geom_tile(color = "white", linewidth = 0.15) +
        scale_fill_manual(values = c("0" = "lightblue", "1" = "firebrick"), na.value = "grey90",
                          breaks = c("1", "0"), labels = c("Methylated", "Unmethylated"), name = "Status") +
        facet_grid(Cluster ~ ., scales = "free_y", space = "free_y", switch = "y") +
        theme_minimal(base_size = 14) +
        theme(axis.text.y = element_blank(), panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, size = 11),
              strip.text.y.left = element_text(angle = 0)) +
        labs(title = paste0(sample_id, ": ASM profile (k-means clusters)"),
             x = "CpG position (sequential)", y = "Reads")
      ggsave(file.path(plot_dir, paste0(safe_name, "_asm.pdf")), p_asm,
             width = 9, height = max(6, min(12, 5 + 0.015 * length(ids))), device = pdf_device)
    }
  }
}

if ("samples" %in% plots) {
  summary_data <- as_data_frame(bundle$summary)
  if (nrow(summary_data) && "Overall_Methylation" %in% names(summary_data)) {
    summary_data <- summary_data %>% filter(is.finite(Overall_Methylation))
    p_batch <- ggplot(summary_data, aes(x = Sample, y = Overall_Methylation, fill = Sample)) +
      geom_col(colour = "black", show.legend = FALSE) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      theme_minimal(base_size = 15) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "PANDA sample summary", x = NULL, y = "Overall methylation (%)")
    ggsave(file.path(plot_dir, "PANDA_sample_summary.pdf"), p_batch,
           width = 11, height = 7, device = pdf_device)
  }
}

if ("group" %in% plots && !is.null(bundle$group_summary)) {
  group_summary <- as_data_frame(bundle$group_summary)
  if (nrow(group_summary)) {
    metric_cols <- c("Mean_Overall_Methylation", "Mean_Amplicon_PDR",
                     "Mean_Window_Epipolymorphism", "Mean_Amplicon_qFDRP")
    sd_cols <- sub("^Mean_", "SD_", metric_cols)
    group_data <- group_summary %>%
      select(Group, N, all_of(metric_cols), all_of(sd_cols))
    mean_data <- group_data %>%
      tidyr::pivot_longer(cols = all_of(metric_cols), names_to = "Metric", values_to = "Mean") %>%
      mutate(Metric = sub("^Mean_", "", Metric))
    sd_data <- group_data %>%
      tidyr::pivot_longer(cols = all_of(sd_cols), names_to = "SD_Metric", values_to = "SD") %>%
      mutate(Metric = sub("^SD_", "", SD_Metric)) %>%
      select(Group, Metric, SD)
    group_data <- mean_data %>%
      left_join(sd_data, by = c("Group", "Metric")) %>%
      mutate(Group_Label = paste0(Group, "\n(n=", N, ")"),
             Metric = recode(Metric,
                             Overall_Overall_Methylation = "Overall methylation (%)",
                             Overall_Methylation = "Overall methylation (%)",
                             Amplicon_PDR = "Amplicon PDR (%)",
                             Window_Epipolymorphism = "Window epipolymorphism",
                             Amplicon_qFDRP = "Amplicon qFDRP"))
    metric_upper <- c("Overall methylation (%)" = 100,
                      "Amplicon PDR (%)" = 100,
                      "Window epipolymorphism" = 1,
                      "Amplicon qFDRP" = 1)
    limit_data <- tidyr::expand_grid(
      Metric = names(metric_upper),
      Group_Label = unique(group_data$Group_Label)
    ) %>%
      mutate(Mean = unname(metric_upper[Metric]), SD = 0, is_limit = TRUE)
    group_data <- group_data %>% mutate(is_limit = FALSE) %>%
      bind_rows(limit_data)
    group_levels <- unique(as.character(group_data$Group))
    palette_values <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
    group_cols <- setNames(rep(palette_values, length.out = length(group_levels)), group_levels)
    p <- ggplot(group_data, aes(x = Group_Label, y = Mean, fill = Group)) +
      geom_blank(data = filter(group_data, is_limit)) +
      geom_col(data = filter(group_data, !is_limit), width = 0.68) +
      geom_errorbar(data = filter(group_data, !is_limit),
                    aes(ymin = Mean - SD, ymax = Mean + SD),
                    width = 0.14, linewidth = 0.5) +
      geom_text(data = filter(group_data, !is_limit),
                aes(label = formatC(Mean, format = "fg", digits = 3)),
                vjust = -0.45, size = 4) +
      facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = group_cols, drop = FALSE) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 11),
            axis.title = element_text(size = 13),
            strip.text = element_text(size = 13),
            panel.grid.minor = element_blank(),
            panel.spacing = grid::unit(2.2, "lines"),
            plot.margin = margin(18, 18, 18, 18)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(title = "PANDA group summary",
           subtitle = "Bars show group means; error bars show SD across samples",
           x = NULL, y = "Mean")
    ggsave(file.path(plot_dir, "PANDA_group_summary.pdf"),
           p, width = 11.5, height = 9.2, device = pdf_device)
    
    ## Also export one clearly labelled PDF per metric for manuscript panels
    ## and presentations where the combined 2x2 layout is too small.
    for (metric_name in names(metric_upper)) {
      one_metric <- group_data %>%
        filter(Metric == metric_name, !is_limit)
      one_plot <- ggplot(one_metric, aes(x = Group_Label, y = Mean, fill = Group)) +
        geom_col(width = 0.62) +
        geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                      width = 0.14, linewidth = 0.5) +
        geom_text(aes(label = formatC(Mean, format = "fg", digits = 3)),
                  vjust = -0.45, size = 4) +
        scale_fill_manual(values = group_cols, drop = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
        coord_cartesian(ylim = c(0, unname(metric_upper[metric_name])),
                        clip = "off") +
        theme_minimal(base_size = 15) +
        theme(legend.position = "none",
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.margin = margin(12, 12, 18, 12)) +
        labs(title = metric_name, x = NULL, y = "Group mean")
      safe_metric <- gsub("[^A-Za-z0-9]+", "_", metric_name)
      safe_metric <- gsub("^_+|_+$", "", safe_metric)
      ggsave(file.path(plot_dir, paste0("PANDA_group_", safe_metric, ".pdf")),
             one_plot, width = 6.5, height = 5.5, device = pdf_device)
    }
  }
}

cat("PANDA plotting completed.\nPlots written to: ", normalizePath(plot_dir), "\n", sep = "")
