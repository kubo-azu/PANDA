#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)


usage <- function() {
  cat(
    paste0(
      "PANDA command-line interface\n\n",
      "Usage:\n",
      "  Rscript cli/panda.R --mode MODE --input INPUT --reference REF [options]\n",
      "  Rscript cli/panda.R --config analysis.json [override options]\n",
      "\n",
      "Required for direct mode:\n",
      "  --mode sanger|amplicon|ngs  --input FILE_OR_DIR  --reference FASTA\n",
      "\n",
      "Common options:\n",
      "  --reference-target NAME  --output-dir DIR  --min-identity N\n",
      "  --min-conversion N  --min-count N  --workers N\n",
      "  --read-mode merged|unmerged  --max-reads N  --max-unique-reads N\n",
      "  --motif-file motifs.txt (one motif per line; all motifs required)\n",
      "  --motifs MOTIF[,MOTIF...]  --ab1-trim-start N  --ab1-trim-end N\n",
      "  --cluster-method kmeans  --k N  --seed N\n",
      "\n",
      "Config files are optional and useful for exact reproducibility.\n",
      "NGS paired-end: --read-mode unmerged (R1/R2 names must contain _R1/_R2 or _1/_2)\n",
      "  Rscript cli/panda.R --help\n"
    )
  )
}


if (!length(args) || "--help" %in% args || "-h" %in% args) {
  usage()
  quit(status = 0L)
}


plot_requested <- "--plot" %in% args
if (plot_requested) {
  stop(
    "Plotting is a separate command. Use: Rscript cli/panda_plot.R --results <output_dir>",
    call. = FALSE
  )
}
plot_requested <- FALSE


arg_value <- function(flag, aliases = character()) {
  idx <- which(args %in% c(flag, aliases))
  if (!length(idx)) return(NULL)
  i <- idx[[1L]]
  if (i == length(args) || startsWith(args[[i + 1L]], "-")) {
    stop(flag, " requires a value.", call. = FALSE)
  }
  args[[i + 1L]]
}

config_path <- arg_value("--config", c("--config-file", "-c"))

script_arg <- grep(
  "^--file=",
  commandArgs(),
  value = TRUE
)

script_path <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  "cli/panda.R"
}

project_root <- normalizePath(
  file.path(dirname(script_path), ".."),
  mustWork = TRUE
)

if (!is.null(config_path)) {
  config_path <- normalizePath(config_path, mustWork = TRUE)
}

mode_override <- NULL
if ("--mode" %in% args) {
  mode_index <- match("--mode", args)
  if (mode_index == length(args)) stop("--mode requires a value.", call. = FALSE)
  mode_override <- tolower(args[[mode_index + 1L]])
}

input_direct <- arg_value("--input", "-i")
reference_direct <- arg_value("--reference", "-r")
reference_target_direct <- arg_value("--reference-target")
min_identity_direct <- arg_value("--min-identity")
min_conversion_direct <- arg_value("--min-conversion")
min_count_direct <- arg_value("--min-count")
workers_direct <- arg_value("--workers")
read_mode_direct <- arg_value("--read-mode")
max_reads_direct <- arg_value("--max-reads")
max_unique_direct <- arg_value("--max-unique-reads")
motifs_direct <- arg_value("--motifs")
ab1_start_direct <- arg_value("--ab1-trim-start")
ab1_end_direct <- arg_value("--ab1-trim-end")
cluster_method_direct <- arg_value("--cluster-method")
k_direct <- arg_value("--k")
seed_direct <- arg_value("--seed")

output_override <- NULL
if ("--output-dir" %in% args) {
  output_index <- match("--output-dir", args)
  if (output_index == length(args)) stop("--output-dir requires a path.", call. = FALSE)
  output_override <- args[[output_index + 1L]]
}

motif_file_override <- NULL
if ("--motif-file" %in% args) {
  motif_file_index <- match("--motif-file", args)
  if (motif_file_index == length(args)) {
    stop("--motif-file requires a path.", call. = FALSE)
  }
  motif_file_override <- args[[motif_file_index + 1L]]
}

core_package <- file.path(project_root, "PANDAcore")

if (!dir.exists(core_package)) {
  stop("PANDAcore directory was not found.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(PANDAcore)
  library(Biostrings)
  library(dplyr)
  library(stringr)
  library(jsonlite)
})


if (!is.null(config_path)) {
  config <- PANDAcore::panda_read_config(config_path)
} else {
  if (is.null(input_direct) || is.null(reference_direct)) {
    usage()
    stop("Without --config, --input and --reference are required.", call. = FALSE)
  }
  config <- list(
    mode = list(if (is.null(mode_override)) "amplicon" else mode_override),
    input = list(input_direct),
    reference = list(reference_direct),
    output = list("panda_cli_results")
  )
  direct_fields <- list(
    reference_target = reference_target_direct,
    min_identity = min_identity_direct,
    min_conversion = min_conversion_direct,
    min_count = min_count_direct,
    workers = workers_direct,
    read_mode = read_mode_direct,
    max_reads = max_reads_direct,
    max_unique_reads = max_unique_direct,
    motifs = motifs_direct,
    ab1_trim_start = ab1_start_direct,
    ab1_trim_end = ab1_end_direct,
    cluster_method = cluster_method_direct,
    k = k_direct,
    seed = seed_direct
  )
  for (nm in names(direct_fields)) {
    if (!is.null(direct_fields[[nm]])) config[[nm]] <- list(direct_fields[[nm]])
  }
}

if (is.null(config$output)) {
  config$output <- list("panda_cli_results")
}
if (!is.null(output_override)) {
  config$output <- list(output_override)
}
if (!is.null(mode_override)) {
  config_mode <- if (!is.null(config$mode)) tolower(as.character(config$mode[[1L]])) else NULL
  if (!is.null(config_mode) && !identical(config_mode, mode_override)) {
    stop(
      "Command-line mode '", mode_override,
      "' conflicts with config mode '", config_mode, "'.",
      call. = FALSE
    )
  }
  config$mode <- list(mode_override)
}
if (!is.null(motif_file_override)) {
  config$motif_file <- list(motif_file_override)
}
PANDAcore::panda_validate_config(config)


mode <- tolower(as.character(config$mode[[1L]]))

if (!mode %in% c("sanger", "amplicon", "ngs")) {
  stop(
    "CLI mode must be one of: sanger, amplicon, ngs.",
    call. = FALSE
  )
}


if (is.null(config$reference)) {
  stop(
    "The configuration must contain a 'reference' field.",
    call. = FALSE
  )
}

reference_path <- as.character(config$reference[[1L]])

if (!file.exists(reference_path)) {
  stop(
    "Reference file does not exist: ",
    reference_path,
    call. = FALSE
  )
}

reference <- Biostrings::readDNAStringSet(reference_path)

if (length(reference) == 0L) {
  stop("The reference FASTA file contains no sequences.", call. = FALSE)
}

reference_target <- if (!is.null(config$reference_target)) {
  as.character(config$reference_target[[1L]])
} else if (length(reference) == 1L) {
  names(reference)[[1L]]
} else {
  stop(
    "The reference FASTA contains multiple targets. Specify 'reference_target' in the configuration.",
    call. = FALSE
  )
}
if (is.na(reference_target) || !nzchar(reference_target) ||
    !(reference_target %in% names(reference))) {
  stop(
    "reference_target was not found in the reference FASTA: ",
    reference_target,
    call. = FALSE
  )
}
reference_sequence <- reference[[reference_target]]


input_files <- PANDAcore::panda_input_files(
  as.character(config$input[[1L]])
)

read_mode <- if (!is.null(config$read_mode)) {
  tolower(as.character(config$read_mode[[1L]]))
} else {
  "merged"
}
if (!read_mode %in% c("merged", "unmerged")) {
  stop("read_mode must be 'merged' or 'unmerged'.", call. = FALSE)
}

input_records <- data.frame(
  sample_id = sub("\\.(fastq|fq|fasta|fa)(\\.gz)?$", "", basename(input_files), ignore.case = TRUE),
  path_r1 = input_files,
  path_r2 = NA_character_,
  stringsAsFactors = FALSE
)
if (identical(mode, "ngs") && identical(read_mode, "unmerged")) {
  base_name <- function(x) sub(
    "(_R[12]|_[12])(_[^.]+)?\\.(fastq|fq|fasta|fa)(\\.gz)?$",
    "", basename(x), ignore.case = TRUE
  )
  direction <- function(x) {
    if (grepl("_R1|_1(_|\\.)", basename(x), ignore.case = TRUE)) return("R1")
    if (grepl("_R2|_2(_|\\.)", basename(x), ignore.case = TRUE)) return("R2")
    NA_character_
  }
  dirs <- vapply(input_files, direction, character(1))
  bases <- vapply(input_files, base_name, character(1))
  paired <- lapply(unique(bases[!is.na(dirs)]), function(bb) {
    r1 <- input_files[bases == bb & dirs == "R1"]
    r2 <- input_files[bases == bb & dirs == "R2"]
    if (!length(r1) || !length(r2)) return(NULL)
    data.frame(sample_id = bb, path_r1 = r1[[1L]], path_r2 = r2[[1L]], stringsAsFactors = FALSE)
  })
  input_records <- bind_rows(paired)
  if (!nrow(input_records)) {
    stop("No paired R1/R2 files were found for read_mode='unmerged'.", call. = FALSE)
  }
}

output_dir <- as.character(config$output[[1L]])

if (!grepl("^(/|[A-Za-z]:[\\\\/])", output_dir)) {
  output_dir <- file.path(project_root, output_dir)
}
output_dir <- normalizePath(output_dir, mustWork = FALSE)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir, call. = FALSE)
}

plot_top_n <- if (!is.null(config$plot_top_n)) {
  as.integer(config$plot_top_n[[1L]])
} else {
  30L
}
if (is.na(plot_top_n) || plot_top_n < 1L) plot_top_n <- 30L
plot_dir <- file.path(output_dir, "plots")
if (plot_requested) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_sample_plots <- function(result, sample_id, plot_dir, top_n = 30L) {
  long_data <- result$long_data
  if (is.null(long_data) || !nrow(long_data)) return(invisible(NULL))
  
  read_counts <- long_data %>%
    distinct(ReadID, Count) %>%
    mutate(Count = ifelse(is.na(Count) | Count < 1, 1L, Count)) %>%
    arrange(desc(Count))
  keep_ids <- head(read_counts$ReadID, top_n)
  plot_data <- long_data %>% filter(ReadID %in% keep_ids)
  plot_counts <- read_counts %>% filter(ReadID %in% keep_ids)
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
  
  read_stats <- plot_data %>%
    group_by(ReadID) %>%
    summarise(
      Mean_Meth = mean(Methylation, na.rm = TRUE) * 100,
      Count = first(Count),
      .groups = "drop"
    )
  
  p_hist <- ggplot(read_stats, aes(x = Mean_Meth, weight = Count)) +
    geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0(sample_id, ": methylation distribution"),
      subtitle = paste0("Top ", nrow(read_stats), " variants for display; metrics use the configured set"),
      x = "Mean methylation per variant (%)", y = "Weighted read count"
    )
  ggsave(file.path(plot_dir, paste0(safe_name, "_methylation_distribution.pdf")),
         p_hist, width = 8, height = 6, device = "pdf")
  
  plot_counts <- plot_counts %>%
    arrange(desc(Count)) %>%
    mutate(ymax = cumsum(Count), ymin = lag(ymax, default = 0))
  plot_data <- plot_data %>% left_join(plot_counts, by = "ReadID")
  p_heat <- ggplot(plot_data) +
    geom_rect(aes(xmin = Position - 2.5, xmax = Position + 2.5,
                  ymin = ymin, ymax = ymax, fill = factor(Methylation)),
              color = NA) +
    scale_fill_manual(values = c("0" = "lightblue", "1" = "firebrick"),
                      labels = c("Unmethylated", "Methylated"), name = "Status") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank()) +
    labs(title = paste0(sample_id, ": abundance heatmap (display Top-N)"),
         x = "CpG position (bp)", y = "Cumulative read count")
  ggsave(file.path(plot_dir, paste0(safe_name, "_abundance_heatmap.pdf")),
         p_heat, width = 9, height = 7, device = "pdf")
  
  plot_data <- plot_data %>%
    mutate(CpG_Index = match(Position, sort(unique(Position))))
  plot_data$ReadID <- factor(plot_data$ReadID, levels = plot_counts$ReadID)
  p_lollipop <- ggplot(plot_data, aes(x = CpG_Index, y = ReadID)) +
    geom_line(aes(group = ReadID), color = "grey80") +
    geom_point(aes(fill = factor(Methylation)), shape = 21, size = 3, color = "black") +
    scale_fill_manual(values = c("0" = "white", "1" = "black"), guide = "none") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = paste0(sample_id, ": methylation lollipop (display Top-N)"),
         x = "CpG index", y = "Variant")
  ggsave(file.path(plot_dir, paste0(safe_name, "_lollipop.pdf")),
         p_lollipop, width = 9, height = 7, device = "pdf")
  invisible(NULL)
}


min_identity <- if (!is.null(config$min_identity)) {
  as.numeric(config$min_identity[[1L]])
} else {
  90
}

min_conversion <- if (!is.null(config$min_conversion)) {
  as.numeric(config$min_conversion[[1L]])
} else {
  95
}

max_reads <- if (!is.null(config$max_reads)) {
  as.integer(config$max_reads[[1L]])
} else {
  Inf
}

max_unique_reads <- if (!is.null(config$max_unique_reads)) {
  as.integer(config$max_unique_reads[[1L]])
} else {
  Inf
}

min_count <- if (!is.null(config$min_count)) {
  as.integer(config$min_count[[1L]])
} else {
  1L
}

ab1_trim_start <- if (!is.null(config$ab1_trim_start)) {
  as.integer(config$ab1_trim_start[[1L]])
} else {
  20L
}
ab1_trim_end <- if (!is.null(config$ab1_trim_end)) {
  as.integer(config$ab1_trim_end[[1L]])
} else {
  20L
}

# Optional in-silico haplotype filter.  A read must contain every motif.
# Motifs are matched literally (not as regular expressions), which is the
# intended behavior for DNA sequence strings and keeps CLI behavior aligned
# with the GUI filter.
motifs <- character()
if (!is.null(config$motifs)) {
  motif_values <- unlist(config$motifs, use.names = FALSE)
  if (length(motif_values) == 1L && grepl(",", motif_values, fixed = TRUE)) {
    motif_values <- strsplit(motif_values, ",", fixed = TRUE)[[1L]]
  }
  motifs <- trimws(as.character(motif_values))
  motifs <- motifs[nzchar(motifs)]
  if (length(motifs) > 0L) {
    motifs <- toupper(motifs)
    if (any(!grepl("^[ACGTNRYKMSWBDHV]+$", motifs))) {
      stop(
        "motifs must contain only IUPAC DNA symbols (A,C,G,T,N,R,Y,K,M,S,W,B,D,H,V).",
        call. = FALSE
      )
    }
  }
}

motif_file <- NULL
if (!is.null(config$motif_file)) {
  motif_file <- as.character(config$motif_file[[1L]])
  if (!grepl("^(/|[A-Za-z]:[\\\\/])", motif_file)) {
    motif_file <- file.path(project_root, motif_file)
  }
  motif_file <- normalizePath(motif_file, mustWork = TRUE)
  file_motifs <- readLines(motif_file, warn = FALSE, encoding = "UTF-8")
  if (length(file_motifs)) {
    file_motifs[[1L]] <- sub("^\\ufeff", "", file_motifs[[1L]])
  }
  file_motifs <- trimws(file_motifs)
  file_motifs <- file_motifs[nzchar(file_motifs) & !startsWith(file_motifs, "#")]
  motifs <- c(motifs, file_motifs)
  motifs <- unique(motifs)
  motifs <- toupper(motifs)
  if (length(motifs) > 0L && any(!grepl("^[ACGTNRYKMSWBDHV]+$", motifs))) {
    stop(
      "motif_file contains a non-DNA motif. Use IUPAC DNA symbols only.",
      call. = FALSE
    )
  }
}

workers <- if (!is.null(config$workers)) {
  as.integer(config$workers[[1L]])
} else {
  16L
}
if (is.na(workers) || workers < 1L) workers <- 1L
workers <- min(workers, 16L)

cluster_method <- if (!is.null(config$cluster_method)) as.character(config$cluster_method[[1L]]) else "kmeans"
cluster_k <- if (!is.null(config$k)) as.integer(config$k[[1L]]) else 2L
cluster_seed <- if (!is.null(config$seed)) as.integer(config$seed[[1L]]) else 11L


read_format <- function(path) {
  lower_path <- tolower(path)
  lower_path <- sub("\\.gz$", "", lower_path)
  
  if (grepl("\\.(fastq|fq)$", lower_path)) {
    return("fastq")
  }
  
  if (grepl("\\.(fasta|fa)$", lower_path)) {
    return("fasta")
  }

  if (grepl("\\.(ab1|abi)$", lower_path)) {
    return("ab1")
  }
  
  stop(
    "Only AB1, FASTA, and FASTQ files are supported: ",
    path,
    call. = FALSE
  )
}


sample_name <- function(path) {
  x <- basename(path)
  sub(
    "\\.(fastq|fq|fasta|fa)(\\.gz)?$",
    "",
    x,
    ignore.case = TRUE
  )
}


summary_rows <- list()
analysis_records <- list()


for (record_index in seq_len(nrow(input_records))) {
  input_file <- input_records$path_r1[[record_index]]
  input_file_r2 <- input_records$path_r2[[record_index]]
  sample_id <- input_records$sample_id[[record_index]]
  format <- read_format(input_file)
  
  cat("Processing: ", sample_id, "\n", sep = "")
  
  reads <- if (identical(format, "ab1")) {
    PANDAcore::process_ab1_files(
      file_paths = input_file,
      file_names = basename(input_file),
      trim_start = ab1_trim_start,
      trim_end = ab1_trim_end
    )
  } else {
    Biostrings::readDNAStringSet(input_file, format = format)
  }
  
  paired_mode <- identical(mode, "ngs") && !is.na(input_file_r2)
  if (paired_mode) {
    reads_r2 <- Biostrings::readDNAStringSet(input_file_r2, format = read_format(input_file_r2))
    n_pair <- min(length(reads), length(reads_r2))
    if (n_pair < 1L) {
      warning("No paired reads found in: ", sample_id)
      next
    }
    reads <- reads[seq_len(n_pair)]
    reads_r2 <- reads_r2[seq_len(n_pair)]
    reads_r2_rc <- as.character(Biostrings::reverseComplement(reads_r2))
    pair_strings <- paste(as.character(reads), reads_r2_rc, sep = "---PAIR---")
    raw_input_n <- length(pair_strings)
    if (length(motifs) > 0L) {
      keep <- vapply(pair_strings, function(x) {
        all(vapply(motifs, function(m) grepl(m, toupper(x), fixed = TRUE), logical(1)))
      }, logical(1))
      pair_strings <- pair_strings[keep]
    }
    motif_retained_n <- length(pair_strings)
    if (!length(pair_strings)) {
      warning("No paired reads remained after filtering: ", sample_id)
      next
    }
    if (is.finite(max_reads) && length(pair_strings) > max_reads) pair_strings <- pair_strings[seq_len(max_reads)]
    pair_counts <- sort(table(pair_strings), decreasing = TRUE)
    pair_counts <- pair_counts[pair_counts >= min_count]
    if (!length(pair_counts)) {
      warning("No paired reads met min_count for: ", sample_id)
      next
    }
    expanded_sequences <- character(2L * length(pair_counts))
    expanded_names <- character(2L * length(pair_counts))
    for (jj in seq_along(pair_counts)) {
      parts <- strsplit(names(pair_counts)[[jj]], "---PAIR---", fixed = TRUE)[[1L]]
      base <- paste0("Rank", jj, "_Count", as.integer(pair_counts[[jj]]))
      ii <- 2L * jj - 1L
      expanded_sequences[c(ii, ii + 1L)] <- parts
      expanded_names[c(ii, ii + 1L)] <- paste0(base, c("_R1", "_R2"))
    }
    reads <- Biostrings::DNAStringSet(expanded_sequences)
    names(reads) <- expanded_names
    raw_read_n <- sum(as.integer(pair_counts))
    dereplicated <- list(reads = reads, counts = rep(as.integer(pair_counts), each = 2L),
                         raw_n = raw_read_n, unique_n = length(pair_counts))
  } else {
    raw_input_n <- length(reads)
    if (length(motifs) > 0L) {
      read_strings <- toupper(as.character(reads))
      keep <- vapply(read_strings, function(x) {
        all(vapply(motifs, function(m) grepl(m, x, fixed = TRUE), logical(1)))
      }, logical(1))
      reads <- reads[keep]
    }
    motif_retained_n <- length(reads)
    if (length(reads) == 0L) {
      warning("No reads found in: ", input_file)
      next
    }
    if (is.finite(max_reads) && length(reads) > max_reads) reads <- reads[seq_len(max_reads)]
    if (is.null(names(reads)) || any(!nzchar(names(reads)))) {
      names(reads) <- paste0(sample_id, "_Read", seq_along(reads))
    }
    raw_read_n <- length(reads)
    if (identical(mode, "sanger")) {
    dereplicated <- list(
      reads = reads,
      counts = rep.int(1L, length(reads)),
      raw_n = raw_read_n,
      unique_n = length(reads)
    )
    } else {
      dereplicated <- PANDAcore::panda_dereplicate_reads(
        reads, min_count = min_count, max_unique_reads = max_unique_reads
      )
      reads <- dereplicated$reads
    }
  }

  if (length(motifs) > 0L) {
    cat(
      "  Motif-filtered reads: ", motif_retained_n,
      if (is.finite(max_reads) && motif_retained_n > raw_read_n) {
        paste0(" (analysis capped at ", raw_read_n, ")")
      } else "",
      "\n", sep = ""
    )
  } else {
    cat(
      "  Input reads: ", raw_input_n,
      if (is.finite(max_reads) && raw_input_n > raw_read_n) {
        paste0(" (analysis capped at ", raw_read_n, ")")
      } else "",
      "\n", sep = ""
    )
  }
  
  cat(
    "  Raw reads: ", raw_read_n,
    "\n",
    "  Unique sequences: ", dereplicated$unique_n,
    "\n",
    sep = ""
  )
  
  result <- PANDAcore::run_bisulfite_alignment(
    genome_seq = reference_sequence,
    reads_set = reads,
    min_identity = min_identity,
    min_conversion = min_conversion,
    return_alignments = TRUE,
    workers = workers
  )
  
  if (is.null(result)) {
    warning("No result was generated for: ", sample_id)
    next
  }
  
  # Alignment is performed on one representative per unique sequence.
  # For unmerged pairs, reproduce GUI behavior by collapsing R1/R2 rows back
  # to one molecule after independent strand-aware alignments.
  if (paired_mode) {
    result$read_summary <- result$read_summary %>%
      mutate(BaseID = sub("_R[12]$", "", ReadID)) %>%
      group_by(BaseID) %>%
      summarise(
        Strand = paste(unique(Strand), collapse = "/"),
        Mismatches = sum(Mismatches), Gaps = sum(Gaps),
        Identity_Pct = mean(Identity_Pct), Meth_Pct = mean(Meth_Pct, na.rm = TRUE),
        Conv_Pct = mean(Conv_Pct), CpG_Count = sum(CpG_Count),
        Pattern = if (any(grepl("excluded", Pattern))) "excluded (Pair Failed QC)" else "Passed",
        .groups = "drop"
      ) %>% rename(ReadID = BaseID)
    result$long_data <- result$long_data %>% mutate(ReadID = sub("_R[12]$", "", ReadID))
    excluded_ids <- result$read_summary$ReadID[grepl("excluded", result$read_summary$Pattern)]
    result$long_data <- result$long_data %>% filter(!ReadID %in% excluded_ids)
    if (!is.null(result$counts$total)) result$counts$total <- result$counts$total / 2
  }

  # Attach the multiplicity so downstream summaries can use read-level weight.
  if (!is.null(result$long_data) && nrow(result$long_data) > 0L) {
    count_ids <- names(reads)
    count_values <- as.integer(dereplicated$counts)
    if (paired_mode) {
      count_ids <- sub("_R[12]$", "", count_ids)
      keep <- !duplicated(count_ids)
      count_ids <- count_ids[keep]
      count_values <- count_values[keep]
    }
    count_map <- data.frame(
      ReadID = count_ids,
      Count = count_values,
      stringsAsFactors = FALSE
    )
    if ("Count" %in% names(result$long_data)) {
      result$long_data$Count <- NULL
    }
    result$long_data <- result$long_data %>%
      left_join(count_map, by = "ReadID")
  }
  
  heterogeneity <- PANDAcore::calculate_heterogeneity(
    result, cluster_method = cluster_method, k = cluster_k, seed = cluster_seed
  )
  
  statistics <- PANDAcore::calculate_quma_stats(
    result,
    mode = if (identical(mode, "sanger")) "Sanger" else "NGS"
  )
  
  if (plot_requested) {
    save_sample_plots(result, sample_id, plot_dir, top_n = plot_top_n)
  }
  
  analysis_records[[sample_id]] <- list(
    Sample = sample_id,
    Input = input_file,
    long_data = result$long_data,
    read_summary = result$read_summary,
    genome_info = result$genome_info,
    alignments = result$alignments,
    heterogeneity = heterogeneity$scores,
    heterogeneity_detail = list(
      meth_mat = data.frame(
        ReadID = rownames(heterogeneity$meth_mat),
        as.data.frame(heterogeneity$meth_mat, check.names = FALSE),
        check.names = FALSE, stringsAsFactors = FALSE
      ),
      clusters = heterogeneity$clusters
    ),
    statistics = statistics$cpg_table
  )
  
  safe_name <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    sample_id
  )
  
  write.csv(
    result$read_summary,
    file = file.path(
      output_dir,
      paste0(safe_name, "_read_summary.csv")
    ),
    row.names = FALSE
  )

  if (!is.null(result$alignments$pairwise)) {
    writeLines(
      unlist(result$alignments$pairwise),
      con = file.path(output_dir, paste0(safe_name, "_alignments.txt"))
    )
  }
  
  write.csv(
    heterogeneity$scores,
    file = file.path(
      output_dir,
      paste0(safe_name, "_heterogeneity.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    statistics$cpg_table,
    file = file.path(
      output_dir,
      paste0(safe_name, "_cpg_stats.csv")
    ),
    row.names = FALSE
  )
  
  passed <- !grepl(
    "^excluded",
    result$read_summary$Pattern
  )
  
  summary_rows[[sample_id]] <- data.frame(
    Sample = sample_id,
    Input = input_file,
    Input_Reads = raw_input_n,
    Reads = raw_read_n,
    Motif_Filtered_Reads = min(motif_retained_n, raw_read_n),
    Motifs = if (length(motifs)) paste(motifs, collapse = ",") else "",
    Motif_File = if (!is.null(motif_file)) motif_file else "",
    Unique_Sequences = dereplicated$unique_n,
    Retained_Reads = sum(dereplicated$counts),
    Min_Count = min_count,
    Aligned_Records = nrow(result$read_summary),
    Passed = sum(passed, na.rm = TRUE),
    Overall_Methylation = statistics$overall,
    Amplicon_PDR = heterogeneity$scores$Value[
      heterogeneity$scores$Metric == "Amplicon PDR"
    ],
    Window_Epipolymorphism = heterogeneity$scores$Value[
      heterogeneity$scores$Metric == "Window Epipolymorphism"
    ],
    Amplicon_qFDRP = heterogeneity$scores$Value[
      heterogeneity$scores$Metric == "Amplicon qFDRP"
    ],
    stringsAsFactors = FALSE
  )
  
  cat(
    "  Reads: ", raw_read_n,
    " (unique: ", dereplicated$unique_n,
    ")\n",
    "  Aligned records: ", nrow(result$read_summary),
    "\n",
    "  Overall methylation: ",
    statistics$overall,
    "\n",
    sep = ""
  )
}


if (length(summary_rows) == 0L) {
  stop("No samples were successfully analyzed.", call. = FALSE)
}

summary_table <- bind_rows(summary_rows)

write.csv(
  summary_table,
  file = file.path(output_dir, "PANDA_summary.csv"),
  row.names = FALSE
)

## Optional replicate-level group summary. Groups are explicit in the
## configuration so that sample naming conventions are never guessed.
group_summary <- data.frame()
if (!is.null(config$groups)) {
  group_rows <- lapply(names(config$groups), function(group_name) {
    members <- as.character(config$groups[[group_name]])
    group_data <- summary_table[summary_table$Sample %in% members, , drop = FALSE]
    if (!nrow(group_data)) {
      warning("No analyzed samples found for group: ", group_name)
      return(NULL)
    }
    metric_names <- c(
      "Overall_Methylation", "Amplicon_PDR",
      "Window_Epipolymorphism", "Amplicon_qFDRP"
    )
    out <- data.frame(
      Group = group_name,
      N = nrow(group_data),
      stringsAsFactors = FALSE
    )
    for (metric in metric_names) {
      values <- group_data[[metric]]
      out[[paste0("Mean_", metric)]] <- mean(values, na.rm = TRUE)
      out[[paste0("SD_", metric)]] <- if (sum(is.finite(values)) > 1L) {
        stats::sd(values, na.rm = TRUE)
      } else {
        NA_real_
      }
    }
    out
  })
  group_summary <- bind_rows(group_rows)
  if (nrow(group_summary)) {
    write.csv(
      group_summary,
      file = file.path(output_dir, "PANDA_group_summary.csv"),
      row.names = FALSE
    )
    if (plot_requested) {
      metric_cols <- c(
        "Mean_Overall_Methylation", "Mean_Amplicon_PDR",
        "Mean_Window_Epipolymorphism", "Mean_Amplicon_qFDRP"
      )
      group_plot_data <- group_summary %>%
        select(Group, all_of(metric_cols)) %>%
        tidyr::pivot_longer(
          cols = all_of(metric_cols),
          names_to = "Metric",
          values_to = "Mean"
        ) %>%
        mutate(Metric = sub("^Mean_", "", Metric))
      p_group <- ggplot(group_plot_data, aes(x = Group, y = Mean, fill = Group)) +
        geom_col() +
        facet_wrap(~ Metric, scales = "free_y") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
        labs(title = "PANDA group means", x = NULL, y = "Mean")
      ggsave(file.path(plot_dir, "PANDA_group_summary.pdf"),
             p_group, width = 10, height = 7, device = "pdf")
    }
  }
}

analysis_bundle <- list(
  schema_version = "0.2",
  generated_at = format(Sys.time(), tz = "UTC"),
  parameters = list(
    mode = mode,
    reference = reference_path,
    reference_target = reference_target,
    read_mode = read_mode,
    min_identity = min_identity,
    min_conversion = min_conversion,
    min_count = min_count,
    motifs = motifs,
    motif_file = motif_file,
    cluster_method = cluster_method,
    k = cluster_k,
    seed = cluster_seed,
    workers = workers
  ),
  summary = summary_table,
  group_summary = group_summary,
  samples = analysis_records
)
jsonlite::write_json(
  analysis_bundle,
  path = file.path(output_dir, "PANDA_analysis.json"),
  dataframe = "rows",
  auto_unbox = TRUE,
  na = "null",
  pretty = TRUE
)

jsonlite::write_json(
  list(
    schema_version = "0.1",
    command = commandArgs(),
    config = config
  ),
  path = file.path(output_dir, "PANDA_run_manifest.json"),
  auto_unbox = TRUE,
  na = "null",
  pretty = TRUE
)

cat(
  "\nPANDA CLI analysis completed.\n",
  "Results written to: ",
  normalizePath(output_dir),
  "\n",
  sep = ""
)
