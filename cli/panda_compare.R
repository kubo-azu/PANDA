#!/usr/bin/env Rscript

## Compare two groups from a completed PANDA CLI analysis.  The calculation is
## delegated to PANDAcore so GUI and CLI use the same CpG-level tests.
args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "PANDA group-comparison command\n\n",
    "Usage:\n",
    "  Rscript cli/panda_compare.R --results results/run\n",
    "  Rscript cli/panda_compare.R --results results/run --group-a sample1,sample2 --group-b sample3,sample4\n",
    "  Rscript cli/panda_compare.R --results results/run --group-a-file groups/wt.txt --group-b-file groups/ko.txt --name-a WT --name-b KO\n\n",
    "Each group file contains one sample ID or input filename per line. Blank lines and # comments are ignored.\n",
    "If group-a/group-b are omitted, groups are read from PANDA_run_manifest.json.\n"
  )
}

if (!length(args) || "--help" %in% args || "-h" %in% args) {
  usage()
  quit(status = 0L)
}

get_arg <- function(flag, required = FALSE) {
  if (!(flag %in% args)) {
    if (required) stop(flag, " is required.", call. = FALSE)
    return(NULL)
  }
  ii <- match(flag, args)
  if (ii == length(args)) stop(flag, " requires a value.", call. = FALSE)
  args[[ii + 1L]]
}

results_dir <- normalizePath(get_arg("--results", required = TRUE), mustWork = TRUE)
json_path <- file.path(results_dir, "PANDA_analysis.json")
manifest_path <- file.path(results_dir, "PANDA_run_manifest.json")
if (!file.exists(json_path)) stop("PANDA_analysis.json was not found in: ", results_dir, call. = FALSE)

suppressPackageStartupMessages({
  library(PANDAcore)
  library(Biostrings)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

bundle <- jsonlite::read_json(json_path, simplifyVector = FALSE)
as_df <- function(x) {
  if (is.null(x) || !length(x)) return(data.frame())
  if (is.data.frame(x)) return(x)
  bind_rows(lapply(x, function(z) as.data.frame(z, stringsAsFactors = FALSE)))
}

sample_records <- bundle$samples
if (is.null(sample_records) || !length(sample_records)) {
  stop("No sample records found in PANDA_analysis.json.", call. = FALSE)
}

reconstruct <- function(x) {
  gi <- x$genome_info
  if (!is.null(gi$cpg_pos)) gi$cpg_pos <- as.numeric(unlist(gi$cpg_pos, use.names = FALSE))
  list(
    long_data = as_df(x$long_data),
    read_summary = as_df(x$read_summary),
    genome_info = gi
  )
}
res_list <- lapply(sample_records, reconstruct)
names(res_list) <- names(sample_records)

manifest <- if (file.exists(manifest_path)) jsonlite::read_json(manifest_path, simplifyVector = FALSE) else NULL
configured_groups <- if (!is.null(manifest$config$groups)) manifest$config$groups else NULL

parse_group <- function(value, label) {
  if (is.null(value)) return(NULL)
  x <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  x <- x[nzchar(x)]
  if (!length(x)) stop(label, " must contain at least one sample.", call. = FALSE)
  x
}

read_group_file <- function(path, label) {
  if (is.null(path)) return(NULL)
  path <- normalizePath(path, mustWork = TRUE)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- trimws(sub("[[:space:]]+#.*$", "", lines))
  lines <- lines[nzchar(lines) & !grepl("^#", lines)]
  if (!length(lines)) stop(label, " is empty: ", path, call. = FALSE)
  lines
}

sample_id_from_entry <- function(x) {
  x <- basename(trimws(x))
  sub("\\.(fastq|fq|fasta|fa|ab1|abi)(\\.gz)?$", "", x, ignore.case = TRUE)
}

resolve_group_members <- function(x, label, available) {
  if (is.null(x)) return(NULL)
  raw <- trimws(as.character(x))
  raw <- raw[nzchar(raw)]
  ids <- vapply(raw, sample_id_from_entry, character(1L))
  out <- unique(c(raw[raw %in% available], ids[ids %in% available]))
  unresolved <- raw[!(raw %in% available | ids %in% available)]
  if (length(unresolved)) {
    stop(label, " contains samples not found in PANDA_analysis.json: ",
         paste(unresolved, collapse = ", "), call. = FALSE)
  }
  if (!length(out)) stop(label, " contains no usable samples.", call. = FALSE)
  out
}

group_a <- parse_group(get_arg("--group-a"), "--group-a")
group_b <- parse_group(get_arg("--group-b"), "--group-b")
group_a_file <- read_group_file(get_arg("--group-a-file"), "--group-a-file")
group_b_file <- read_group_file(get_arg("--group-b-file"), "--group-b-file")
if (xor(is.null(group_a_file), is.null(group_b_file))) {
  stop("--group-a-file and --group-b-file must be supplied together.", call. = FALSE)
}
if (!is.null(group_a_file) || !is.null(group_b_file)) {
  if (!is.null(group_a) || !is.null(group_b)) {
    stop("Use either --group-a/--group-b or --group-a-file/--group-b-file, not both.", call. = FALSE)
  }
  group_a <- group_a_file
  group_b <- group_b_file
}
group_names <- c(get_arg("--name-a"), get_arg("--name-b"))

if (is.null(group_a) || is.null(group_b)) {
  if (is.null(configured_groups) || length(configured_groups) < 2L) {
    stop("Provide --group-a and --group-b, or define at least two groups in the analysis config.", call. = FALSE)
  }
  group_a <- as.character(configured_groups[[1L]])
  group_b <- as.character(configured_groups[[2L]])
  group_names <- names(configured_groups)[1:2]
}
if (length(group_names) != 2L || any(is.null(group_names)) || any(!nzchar(group_names))) {
  group_names <- c("Group 1", "Group 2")
}
group_a <- resolve_group_members(group_a, "Group A", names(res_list))
group_b <- resolve_group_members(group_b, "Group B", names(res_list))
overlap <- intersect(group_a, group_b)
if (length(overlap)) {
  stop("A sample cannot belong to both groups: ", paste(overlap, collapse = ", "), call. = FALSE)
}
missing_samples <- setdiff(c(group_a, group_b), names(res_list))
if (length(missing_samples)) stop("Samples not found in analysis: ", paste(missing_samples, collapse = ", "), call. = FALSE)

reference_path <- bundle$parameters$reference
reference_target <- bundle$parameters$reference_target
if (is.null(reference_path) || is.null(reference_target)) {
  stop("Reference metadata is missing from PANDA_analysis.json.", call. = FALSE)
}
reference <- Biostrings::readDNAStringSet(reference_path)
if (!(reference_target %in% names(reference))) stop("Reference target not found: ", reference_target, call. = FALSE)

comparison <- PANDAcore::analyze_group_comparison(
  res_list[group_a], res_list[group_b], reference[[reference_target]],
  g1_name = group_names[[1L]], g2_name = group_names[[2L]]
)
if (is.null(comparison)) stop("The two groups did not contain usable results.", call. = FALSE)

safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
out_prefix <- file.path(
  results_dir,
  paste0("PANDA_comparison_", safe(group_names[[1L]]), "_vs_", safe(group_names[[2L]]))
)
site_path <- paste0(out_prefix, "_sites.csv")
stats_path <- paste0(out_prefix, "_statistics.json")
diff_path <- paste0(out_prefix, "_difference.pdf")
bar_path <- paste0(out_prefix, "_summary.pdf")

write.csv(comparison$site_table, site_path, row.names = FALSE)
jsonlite::write_json(
  list(
    group_a = group_names[[1L]], group_b = group_names[[2L]],
    samples_a = group_a, samples_b = group_b,
    replicate_p_used = comparison$replicate_p_used,
    overall_test = comparison$overall_test,
    pooled_overall_test = comparison$pooled_overall_test,
    site_test_uses_replicates = comparison$site_test_uses_replicates,
    pooled_u_test_p = comparison$pooled_u_test_p,
    u_test_p = comparison$u_test_p,
    sample_values = comparison$sample_values,
    sample_summary = comparison$sample_summary,
    summary = comparison$summary,
    site_table = comparison$site_table
  ), stats_path, dataframe = "rows", auto_unbox = TRUE, na = "null", pretty = TRUE
)

site <- comparison$site_table
site$Delta <- site$Pct_1 - site$Pct_2
p_diff <- ggplot(site, aes(x = factor(Position), y = Delta, fill = Delta > 0)) +
  geom_col(width = 0.72, na.rm = TRUE) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"), guide = "none") +
  geom_hline(yintercept = 0, colour = "grey35") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = paste0("Methylation difference: ", group_names[[1L]], " vs ", group_names[[2L]]),
       x = "CpG position (bp)", y = "Difference in methylation (%)")
ggsave(diff_path, p_diff, width = 9, height = 6, device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")

bar <- comparison$sample_summary
bar$Group <- factor(bar$Group, levels = group_names)
bar$Label_Y <- ifelse(bar$Mean > 92, bar$Mean - 4, bar$Mean + 4)
sample_values <- comparison$sample_values
sample_values$Group <- factor(sample_values$Group, levels = group_names)
p_bar <- ggplot(bar, aes(x = Group, y = Mean, fill = Group)) +
  geom_col(width = 0.65) +
  geom_errorbar(
    aes(ymin = pmax(0, Mean - SD), ymax = pmin(100, Mean + SD)),
    width = 0.14, linewidth = 0.5, na.rm = TRUE
  ) +
  geom_point(
    data = sample_values,
    aes(x = Group, y = Overall_Methylation),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.07, height = 0, seed = 11),
    shape = 21, size = 3, stroke = 0.7, fill = "white", colour = "black"
  ) +
  scale_fill_manual(values = c("#0072B2", "#D55E00"), drop = FALSE) +
  geom_text(
    aes(y = Label_Y, label = sprintf("%.2f", Mean)),
    vjust = 0.5, size = 4
  ) +
  scale_y_continuous(breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.08))) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(
    title = "Group methylation summary",
    subtitle = "Bars: sample means; error bars: SD; points: samples",
    x = NULL, y = "Mean methylation (%)"
  )
ggsave(bar_path, p_bar, width = 7.2, height = 5.5, device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")

cat("PANDA comparison completed.\n")
cat("Sites: ", site_path, "\nStatistics: ", stats_path, "\n", sep = "")
cat("Plots: ", diff_path, " and ", bar_path, "\n", sep = "")
