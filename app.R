library(shiny)
library(shinyjs)
library(Biostrings)
library(tidyverse)
library(DT)
library(bslib)
library(patchwork)
library(sangerseqR)

# Bioconductor provides pairwise alignment helpers through pwalign.  It is a
# required dependency for the current app; the deprecated Biostrings alignment
# fallback is intentionally not used.
if (!requireNamespace("pwalign", quietly = TRUE)) {
  stop("Package 'pwalign' is required. Install it through the project environment.")
}
if (!requireNamespace("PANDAcore", quietly = TRUE)) {
  stop("Package 'PANDAcore' is required. Install it with renv::install('./PANDAcore').")
}
.panda_has_pwalign <- TRUE
.panda_alignment_cache <- new.env(parent = emptyenv())
.panda_alignment_cache_limit <- 10000L
.panda_pairwiseAlignment <- function(...) {
  args <- list(...)
  pattern_key <- as.character(args$pattern)
  subject_key <- as.character(args$subject)
  type_key <- if (is.null(args$type)) "" else as.character(args$type)
  gap_open_key <- if (is.null(args$gapOpening)) "" else as.character(args$gapOpening)
  gap_ext_key <- if (is.null(args$gapExtension)) "" else as.character(args$gapExtension)
  key <- paste(type_key, gap_open_key, gap_ext_key, pattern_key, subject_key, sep = "\r")
  if (exists(key, envir = .panda_alignment_cache, inherits = FALSE)) {
    return(get(key, envir = .panda_alignment_cache, inherits = FALSE))
  }
  aln <- do.call(pwalign::pairwiseAlignment, args)
  if (length(ls(envir = .panda_alignment_cache)) < .panda_alignment_cache_limit) {
    assign(key, aln, envir = .panda_alignment_cache)
  }
  aln
}
.panda_nucleotideSubstitutionMatrix <- function(...) {
  pwalign::nucleotideSubstitutionMatrix(...)
}

# ==============================================================================
# 0. GLOBAL SETTINGS
# ==============================================================================
# Maximum 300MB for large NGS FASTQ files
options(shiny.maxRequestSize = 300 * 1024^2)
# Shiny reactive contexts cannot safely cross forked processes.  The shared
# engine therefore uses the portable socket backend for GUI requests.
Sys.setenv(PANDA_NO_FORK = "1")
options(panda.no_fork = TRUE)

# Optional NGS parallelism. The alignment algorithm and results are unchanged;
# workers only distribute independent representative alignments. Users can
# override this in the GUI or with PANDA_WORKERS before launching the app.
# The public default and maximum are both 16; constrained hosts can override
# PANDA_WORKERS (for example, use 2 on a free Hugging Face Space).
.panda_default_workers <- 16L
.panda_workers <- suppressWarnings(as.integer(Sys.getenv("PANDA_WORKERS", as.character(.panda_default_workers))))
if (is.na(.panda_workers) || .panda_workers < 1L) .panda_workers <- 1L

# ==============================================================================
# 1. CORE LOGIC (Backend Functions)
# ==============================================================================

run_bisulfite_alignment <- function(genome_seq, reads_set, 
                                    min_identity = 90, 
                                    min_conversion = 95,
                                    return_alignments = TRUE,
                                    workers = 16L) {
  
  if (is(genome_seq, "DNAStringSet")) genome_seq <- genome_seq[[1]]
  genome_seq <- DNAString(as.character(genome_seq))
  
  # ----------------------------------------------------------------------------
  # Alignment Scoring Parameters
  # ----------------------------------------------------------------------------
  sub_mat <- .panda_nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = FALSE)
  gap_op <- -10  # Biostrings default for strict gap opening
  gap_ext <- -4  # Biostrings default for strict gap extension
  
  # ----------------------------------------------------------------------------
  # Reference Auto-Correction (Strand Detection)
  # ----------------------------------------------------------------------------
  was_flipped <- FALSE
  
  # Step 1: Take a small sample of top reads (up to 10) for a quick strand detection test.
  n_test <- min(10, length(reads_set))
  
  if (n_test > 0) {
    test_reads <- reads_set[1:n_test]
    
    test_conv <- function(test_gen) {
      test_gen_char <- as.character(test_gen)
      gen_T <- chartr("C", "T", test_gen_char)
      c_conv <- 0; c_unconv <- 0
      
      for(i in seq_len(n_test)) {
        r <- as.character(test_reads[[i]])
        
        r_F <- r; r_F_T <- chartr("C", "T", r_F)
        aln_F <- .panda_pairwiseAlignment(pattern = r_F_T, subject = gen_T, type = "global", 
                                          substitutionMatrix = sub_mat, gapOpening = gap_op, gapExtension = gap_ext)
        
        r_R <- as.character(reverseComplement(DNAString(r)))
        r_R_T <- chartr("C", "T", r_R)
        aln_R <- .panda_pairwiseAlignment(pattern = r_R_T, subject = gen_T, type = "global", 
                                          substitutionMatrix = sub_mat, gapOpening = gap_op, gapExtension = gap_ext)
        
        if (score(aln_F) >= score(aln_R)) {
          best_aln <- aln_F; final_r <- r_F
        } else {
          best_aln <- aln_R; final_r <- r_R
        }
        
        aln_pat_converted <- as.character(pattern(best_aln))
        raw_bases <- strsplit(final_r, "")[[1]]
        aln_template <- strsplit(aln_pat_converted, "")[[1]]
        
        reconstructed <- character(length(aln_template))
        ## Preserve the original read coordinate when restoring C/T states.
        raw_idx <- start(pattern(best_aln))
        for(j in seq_along(aln_template)) {
          if(aln_template[j] == "-") reconstructed[j] <- "-"
          else {
            if(raw_idx <= length(raw_bases)) { reconstructed[j] <- raw_bases[raw_idx]; raw_idx <- raw_idx + 1 }
            else reconstructed[j] <- "N"
          }
        }
        
        aln_sub_str <- as.character(subject(best_aln))
        start_genome <- start(subject(best_aln))
        curr_g_pos <- start_genome - 1
        sub_chars <- strsplit(aln_sub_str, "")[[1]]
        
        for(k in seq_along(sub_chars)) {
          s_char <- sub_chars[k]; p_char <- reconstructed[k]
          if(s_char != "-") curr_g_pos <- curr_g_pos + 1
          if(s_char == "-" || p_char == "-") next
          
          orig_g_base <- substring(test_gen_char, curr_g_pos, curr_g_pos)
          
          if(orig_g_base == "C") {
            if(p_char == "T") c_conv <- c_conv + 1
            else if(p_char == "C") c_unconv <- c_unconv + 1
          }
        }
      }
      if((c_conv + c_unconv) == 0) return(0) else return((c_conv / (c_conv + c_unconv)) * 100)
    }
    
    conv_fwd <- test_conv(genome_seq)
    conv_rc <- test_conv(reverseComplement(genome_seq))
    
    if(conv_rc > 50 && conv_rc > (conv_fwd + 20)) {
      genome_seq <- reverseComplement(genome_seq)
      was_flipped <- TRUE
    }
  }
  
  # ----------------------------------------------------------------------------
  # Main Alignment and CpG Calling
  # ----------------------------------------------------------------------------
  genome_seq_char <- as.character(genome_seq)
  cpg_hits <- matchPattern("CG", genome_seq_char)
  cpg_sites <- start(cpg_hits) 
  genome_T <- chartr("C", "T", genome_seq_char)
  
  res_ids <- character(); res_strands <- character()
  res_pos <- integer(); res_meth <- integer()
  
  sum_id <- character(); sum_strand <- character()
  sum_mm <- integer(); sum_gaps <- integer()
  sum_ident <- numeric(); sum_meth_pct <- numeric()
  sum_conv_pct <- numeric(); sum_pattern <- character()
  sum_meth_cpgs <- integer()
  
  pairwise_txt_list <- list(); multi_align_seqs <- list()
  
  n_total_reads <- length(reads_set)
  n_excluded <- 0 
  
  workers <- suppressWarnings(as.integer(workers))
  if (is.na(workers) || workers < 1L) workers <- 1L
  if (workers > 1L && .Platform$OS.type == "windows") {
    message("  Windows detected: using serial alignment for portability.")
    workers <- 1L
  }
  
  align_one_read <- function(i) {
    raw_read_seq <- reads_set[[i]]
    read_F <- raw_read_seq
    read_F_T <- chartr("C", "T", read_F)
    aln_F <- .panda_pairwiseAlignment(
      pattern = read_F_T, subject = genome_T, type = "global",
      substitutionMatrix = sub_mat, gapOpening = gap_op,
      gapExtension = gap_ext
    )
    read_R <- reverseComplement(raw_read_seq)
    read_R_T <- chartr("C", "T", read_R)
    aln_R <- .panda_pairwiseAlignment(
      pattern = read_R_T, subject = genome_T, type = "global",
      substitutionMatrix = sub_mat, gapOpening = gap_op,
      gapExtension = gap_ext
    )
    list(read_F = read_F, read_R = read_R, aln_F = aln_F, aln_R = aln_R)
  }
  
  alignment_pairs <- if (workers > 1L && n_total_reads > 1L) {
    message("  Parallel alignment: ", n_total_reads,
            " representatives on ", workers, " workers")
    parallel::mclapply(
      seq_len(n_total_reads), align_one_read,
      mc.cores = workers, mc.preschedule = TRUE
    )
  } else {
    lapply(seq_len(n_total_reads), align_one_read)
  }
  
  for (i in seq_len(n_total_reads)) {
    if (i == 1L || i %% 100L == 0L || i == n_total_reads) {
      message("  Aligning representative ", i, "/", n_total_reads)
    }
    if (i %% 500 == 0) gc()
    
    alignment_pair <- alignment_pairs[[i]]
    raw_read_seq <- alignment_pair$read_F
    read_name <- names(reads_set)[i]
    
    read_F <- alignment_pair$read_F
    read_R <- alignment_pair$read_R
    aln_F <- alignment_pair$aln_F
    aln_R <- alignment_pair$aln_R
    
    if (score(aln_F) >= score(aln_R)) {
      best_aln <- aln_F; final_read_seq <- read_F; strand <- "Forward"
    } else {
      best_aln <- aln_R; final_read_seq <- read_R; strand <- "Reverse"
    }
    
    aln_pat_converted <- as.character(pattern(best_aln))
    raw_bases <- strsplit(as.character(final_read_seq), "")[[1]]
    aln_template <- strsplit(aln_pat_converted, "")[[1]]
    
    reconstructed_pat_chars <- character(length(aln_template))
    ## Preserve the raw-read coordinate when restoring C/T states.
    raw_idx <- start(pattern(best_aln))
    
    for(j in seq_along(aln_template)) {
      if(aln_template[j] == "-") {
        reconstructed_pat_chars[j] <- "-"
      } else {
        if(raw_idx <= length(raw_bases)) {
          reconstructed_pat_chars[j] <- raw_bases[raw_idx]
          raw_idx <- raw_idx + 1
        } else {
          reconstructed_pat_chars[j] <- "N"
        }
      }
    }
    
    aln_sub_str <- as.character(subject(best_aln))
    aln_len <- length(reconstructed_pat_chars)
    n_match <- nmatch(best_aln)
    
    n_gaps_read <- str_count(aln_pat_converted, "-")
    n_gaps_genome <- str_count(aln_sub_str, "-")
    total_gaps <- n_gaps_read + n_gaps_genome
    identity_score <- (n_match / aln_len) * 100
    
    start_genome <- start(subject(best_aln))
    curr_g_pos <- start_genome - 1
    
    sub_chars <- strsplit(aln_sub_str, "")[[1]]
    pat_chars <- reconstructed_pat_chars
    
    tmp_meth <- integer(); tmp_pos <- integer()
    conv_C_count <- 0; unconv_C_count <- 0
    
    for (k in 1:aln_len) {
      s_char <- sub_chars[k]; p_char <- pat_chars[k]
      
      if (s_char != "-") curr_g_pos <- curr_g_pos + 1
      if (s_char == "-" || p_char == "-") next
      
      ## genome_seq_char is a single character string; use substring() for
      ## one-based genomic coordinates (subseq() treats it as length 1).
      orig_g_base <- substring(genome_seq_char, curr_g_pos, curr_g_pos)
      
      if (orig_g_base == "C") {
        if (curr_g_pos %in% cpg_sites) {
          if (p_char == "C") { tmp_meth <- c(tmp_meth, 1L); tmp_pos <- c(tmp_pos, curr_g_pos) }
          else if (p_char == "T") { tmp_meth <- c(tmp_meth, 0L); tmp_pos <- c(tmp_pos, curr_g_pos) }
        } else {
          if (p_char == "C") unconv_C_count <- unconv_C_count + 1
          else if (p_char == "T") conv_C_count <- conv_C_count + 1
        }
      }
    }
    
    total_cph <- unconv_C_count + conv_C_count
    conv_rate <- if (total_cph > 0) (conv_C_count / total_cph) * 100 else 100
    
    exclusion_reason <- ""
    if (identity_score < min_identity) exclusion_reason <- paste0("excluded (Id:", round(identity_score,1), "%)")
    else if (conv_rate < min_conversion) exclusion_reason <- paste0("excluded (Conv:", round(conv_rate,1), "%)")
    
    if (exclusion_reason != "") n_excluded <- n_excluded + 1
    
    sum_id <- c(sum_id, read_name); sum_strand <- c(sum_strand, strand)
    sum_mm <- c(sum_mm, aln_len - n_match); sum_gaps <- c(sum_gaps, total_gaps)
    sum_ident <- c(sum_ident, round(identity_score, 1))
    
    if (exclusion_reason == "") {
      m_pct <- if(length(tmp_meth)>0) round(mean(tmp_meth)*100, 1) else NA
      sum_meth_pct <- c(sum_meth_pct, m_pct)
      sum_conv_pct <- c(sum_conv_pct, round(conv_rate, 1))
      sum_pattern <- c(sum_pattern, "Passed")
      sum_meth_cpgs <- c(sum_meth_cpgs, length(tmp_meth))
      
      if (length(tmp_meth) > 0) {
        n_sites <- length(tmp_meth)
        res_ids <- c(res_ids, rep(read_name, n_sites))
        res_strands <- c(res_strands, rep(strand, n_sites))
        res_pos <- c(res_pos, tmp_pos)
        res_meth <- c(res_meth, tmp_meth)
      }
      
      if (return_alignments) {
        pair_txt <- paste0("> ", read_name, "\nGen: ", aln_sub_str, "\nSeq: ", paste(reconstructed_pat_chars, collapse=""), "\n")
        pairwise_txt_list[[length(pairwise_txt_list)+1]] <- pair_txt
      }
    } else {
      sum_meth_pct <- c(sum_meth_pct, NA)
      sum_conv_pct <- c(sum_conv_pct, round(conv_rate, 1))
      sum_pattern <- c(sum_pattern, exclusion_reason)
      sum_meth_cpgs <- c(sum_meth_cpgs, 0)
    }
  }
  
  if (length(sum_id) == 0) return(NULL)
  
  read_summary_df <- data.frame(
    ReadID = sum_id, Strand = sum_strand, Mismatches = sum_mm, Gaps = sum_gaps,
    Identity_Pct = sum_ident, Meth_Pct = sum_meth_pct, Conv_Pct = sum_conv_pct,
    CpG_Count = sum_meth_cpgs,
    Pattern = sum_pattern, stringsAsFactors = FALSE
  )
  
  long_data_df <- data.frame(
    ReadID = res_ids, Strand = res_strands, Position = res_pos, Methylation = res_meth,
    stringsAsFactors = FALSE
  )
  
  return(list(
    long_data = long_data_df,
    read_summary = read_summary_df,
    genome_info = list(len = nchar(genome_seq_char), n_cpg = length(cpg_sites), cpg_pos = cpg_sites, seq = genome_seq_char),
    counts = list(total = n_total_reads, used = n_total_reads - n_excluded, excluded = n_excluded),
    alignments = list(pairwise = pairwise_txt_list, multi = multi_align_seqs),
    was_flipped = was_flipped
  ))
}

process_ab1_files <- function(file_paths, file_names, trim_start = 20, trim_end = 20) {
  seq_list <- DNAStringSet()
  for (i in seq_along(file_paths)) {
    tryCatch({
      sanger <- readsangerseq(file_paths[i]); seq <- primarySeq(sanger)
      if (length(seq) > (trim_start + trim_end)) seq <- subseq(seq, start=trim_start+1, end=length(seq)-trim_end)
      current_set <- DNAStringSet(seq); names(current_set) <- file_names[i]
      seq_list <- c(seq_list, current_set)
    }, error = function(e) warning(paste("Failed:", file_names[i])))
  }
  return(seq_list)
}

# Use the shared PANDAcore engine for all GUI alignments.  This keeps the GUI
# and CLI on the same OS-independent BiocParallel backend and avoids maintaining
# a second alignment implementation in app.R.
run_bisulfite_alignment <- function(genome_seq, reads_set,
                                    min_identity = 90,
                                    min_conversion = 95,
                                    return_alignments = TRUE,
                                    workers = 16L) {
  ## Force all Shiny inputs in the parent process before creating workers.
  ## Otherwise a lazy promise such as input$ngs_ident can be evaluated inside
  ## a worker, where no reactive context exists.
  genome_seq <- force(genome_seq)
  reads_set <- force(reads_set)
  min_identity <- as.numeric(force(min_identity))
  min_conversion <- as.numeric(force(min_conversion))
  return_alignments <- isTRUE(force(return_alignments))
  workers <- as.integer(force(workers))
  PANDAcore::run_bisulfite_alignment(
    genome_seq = genome_seq,
    reads_set = reads_set,
    min_identity = min_identity,
    min_conversion = min_conversion,
    return_alignments = return_alignments,
    workers = workers
  )
}

## Read one FASTA/FASTQ file while preserving multi-line FASTA records and
## validating the four-line FASTQ structure. Returned sequences are uppercase.
.panda_read_sequences <- function(filepath, filename, format = c("fastq", "fasta")) {
  format <- match.arg(format)
  if (is.na(filepath) || is.na(filename) || !file.exists(filepath)) return(character())
  con <- if (grepl("\\.gz$", filename, ignore.case = TRUE)) gzfile(filepath, "rt") else file(filepath, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  if (!length(lines)) return(character())
  if (format == "fastq") {
    if (length(lines) %% 4L != 0L || any(!startsWith(lines[seq(1L, length(lines), by = 4L)], "@")))
      stop("Invalid FASTQ structure: expected 4 lines per record.")
    seqs <- lines[seq(2L, length(lines), by = 4L)]
    quals <- lines[seq(4L, length(lines), by = 4L)]
    if (any(nchar(seqs) != nchar(quals))) stop("FASTQ sequence/quality length mismatch.")
  } else {
    hdr <- which(startsWith(lines, ">"))
    if (!length(hdr)) stop("FASTA file contains no header.")
    ends <- c(hdr[-1L] - 1L, length(lines))
    seqs <- vapply(seq_along(hdr), function(ii) {
      x <- lines[(hdr[ii] + 1L):ends[ii]]
      paste0(x[!startsWith(x, ">")], collapse = "")
    }, character(1))
  }
  seqs <- toupper(gsub("\\s+", "", seqs))
  if (any(!grepl("^[ACGTN]+$", seqs))) stop("Sequence contains unsupported bases.")
  seqs[nzchar(seqs)]
}

.panda_expand_dereplicated <- function(counts, top_n = Inf, paired = FALSE) {
  if (!length(counts)) return(DNAStringSet())
  counts <- head(sort(counts, decreasing = TRUE), top_n)
  if (!paired) {
    out <- DNAStringSet(names(counts))
    names(out) <- paste0("Rank", seq_along(counts), "_Count", as.integer(counts))
    return(out)
  }
  seqs <- character(); ids <- character()
  for (ii in seq_along(counts)) {
    parts <- strsplit(names(counts)[ii], "---PAIR---", fixed = TRUE)[[1L]]
    base <- paste0("Rank", ii, "_Count", as.integer(counts[ii]))
    seqs <- c(seqs, parts[1L], parts[2L])
    ids <- c(ids, paste0(base, "_R1"), paste0(base, "_R2"))
  }
  out <- DNAStringSet(seqs); names(out) <- ids; out
}

.panda_finalize_ngs_result <- function(res, is_unmerged = FALSE) {
  if (is.null(res)) return(NULL)
  if (is_unmerged) {
    res$read_summary <- res$read_summary %>%
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
    res$long_data <- res$long_data %>% mutate(ReadID = sub("_R[12]$", "", ReadID))
    excluded_ids <- res$read_summary$ReadID[grepl("excluded", res$read_summary$Pattern)]
    res$long_data <- res$long_data %>% filter(!ReadID %in% excluded_ids)
    res$counts$total <- res$counts$total / 2
  }
  res$long_data <- res$long_data %>%
    mutate(Count = as.integer(str_extract(ReadID, "(?<=Count)\\d+")),
           Count = replace_na(Count, 1L))
  res$read_summary <- res$read_summary %>%
    mutate(Count = as.integer(str_extract(ReadID, "(?<=Count)\\d+")),
           Count = replace_na(Count, 1L))
  res
}

.panda_weighted_qfdrp <- function(meth_mat, weights, min_shared_cpg = 1L) {
  n <- nrow(meth_mat)
  if (n < 2L) return(0)
  weights <- as.numeric(weights)
  weights[is.na(weights) | weights < 1] <- 1
  masks <- apply(!is.na(meth_mat), 1L, paste0, collapse = "")
  groups <- split(seq_len(n), masks)
  q_sum <- 0; q_den <- 0
  for (aa in seq_along(groups)) {
    ia <- groups[[aa]]; xa <- meth_mat[ia, , drop = FALSE]; wa <- weights[ia]
    shared_self <- !is.na(meth_mat[ia[[1L]], ]); k_self <- sum(shared_self)
    if (k_self >= min_shared_cpg && length(ia) > 1L) {
      total_a <- sum(wa)
      meth_a <- colSums(xa[, shared_self, drop = FALSE] * wa)
      q_sum <- q_sum + sum(meth_a * (total_a - meth_a)) / k_self
      q_den <- q_den + (total_a^2 - sum(wa^2)) / 2
    }
    if (aa == length(groups)) next
    for (bb in (aa + 1L):length(groups)) {
      ib <- groups[[bb]]; xb <- meth_mat[ib, , drop = FALSE]; wb <- weights[ib]
      shared <- !is.na(meth_mat[ia[[1L]], ]) & !is.na(meth_mat[ib[[1L]], ])
      k_shared <- sum(shared)
      if (k_shared < min_shared_cpg) next
      total_a <- sum(wa); total_b <- sum(wb)
      meth_a <- colSums(xa[, shared, drop = FALSE] * wa)
      meth_b <- colSums(xb[, shared, drop = FALSE] * wb)
      q_sum <- q_sum + sum(meth_a * (total_b - meth_b) +
                             (total_a - meth_a) * meth_b) / k_shared
      q_den <- q_den + total_a * total_b
    }
  }
  if (q_den > 0) q_sum / q_den else 0
}

calculate_heterogeneity <- function(res_obj, min_shared_cpg = 1L,
                                    min_window_coverage = 1L,
                                    cluster_method = "kmeans",
                                    k = 2L,
                                    seed = 11L) {
  df_long <- res_obj$long_data
  cpg_sites <- res_obj$genome_info$cpg_pos
  empty_scores <- data.frame(
    Metric = c("Amplicon PDR", "Window Epipolymorphism", "Amplicon qFDRP"),
    Value = c(0, 0, 0)
  )
  if (is.null(df_long) || nrow(df_long) == 0 || length(cpg_sites) == 0)
    return(list(scores = empty_scores, meth_mat = matrix(NA),
                clusters = NULL, pdr_by_cpg = data.frame(),
                epipolymorphism_by_window = data.frame(),
                qfdrp_by_cpg = data.frame()))
  
  ## One row per dereplicated molecule and an explicit abundance column.
  ## Sanger rows have Count=1; NGS rows retain their dereplication count.
  if (!"Count" %in% names(df_long)) df_long$Count <- 1L
  df_long$Count[is.na(df_long$Count) | df_long$Count < 1] <- 1L
  molecule_counts <- df_long %>%
    distinct(ReadID, Count) %>%
    group_by(ReadID) %>% summarise(Count = max(Count), .groups = "drop")
  meth_mat <- df_long %>%
    select(ReadID, Position, Methylation) %>%
    distinct(ReadID, Position, .keep_all = TRUE) %>%
    pivot_wider(names_from = Position, values_from = Methylation) %>%
    column_to_rownames("ReadID")
  miss <- setdiff(as.character(cpg_sites), colnames(meth_mat))
  if (length(miss) > 0) for (cc in miss) meth_mat[[cc]] <- NA_real_
  meth_mat <- meth_mat[, as.character(cpg_sites), drop = FALSE]
  molecule_counts <- molecule_counts[match(rownames(meth_mat), molecule_counts$ReadID), ]
  molecule_counts$Count[is.na(molecule_counts$Count)] <- 1L
  weights <- molecule_counts$Count
  
  ## Amplicon-level PDR: the same discordance concept as the reference,
  ## aggregated over molecules that contain at least four observed CpGs.
  observed_n <- rowSums(!is.na(meth_mat))
  discordant <- apply(meth_mat, 1, function(x) {
    x <- x[!is.na(x)]
    length(x) >= 4L && length(unique(x)) > 1L
  })
  eligible <- observed_n >= 4L
  pdr_score <- if (any(eligible))
    100 * sum(weights[eligible] * discordant[eligible]) / sum(weights[eligible]) else 0
  
  ## Epipolymorphism is calculated independently for each consecutive
  ## four-CpG window, then summarized across windows. Counts are weighted.
  epi_rows <- list()
  if (length(cpg_sites) >= 4L) {
    for (ww in seq_len(length(cpg_sites) - 3L)) {
      x <- meth_mat[, ww:(ww + 3L), drop = FALSE]
      keep <- complete.cases(x)
      if (sum(weights[keep]) < min_window_coverage || !any(keep)) next
      patterns <- apply(x[keep, , drop = FALSE], 1, paste0, collapse = "")
      tab <- tapply(weights[keep], patterns, sum)
      pk <- tab / sum(tab)
      epi_rows[[length(epi_rows) + 1L]] <- data.frame(
        Window = ww, Start_CpG = cpg_sites[ww], End_CpG = cpg_sites[ww + 3L],
        Coverage = sum(weights[keep]), Epipolymorphism = 1 - sum(pk^2)
      )
    }
  }
  epi_df <- if (length(epi_rows)) bind_rows(epi_rows) else data.frame()
  epipoly_score <- if (nrow(epi_df)) mean(epi_df$Epipolymorphism) else 0
  
  ## Amplicon qFDRP: abundance-weighted normalized Hamming distances
  ## over the CpGs shared by each pair of dereplicated molecules. This is
  ## intentionally region-level for long-range amplicon phasing.
  qfdrp_score <- .panda_weighted_qfdrp(
    meth_mat, weights, min_shared_cpg = min_shared_cpg
  )

  ## Clustering is used only to order the ASM heatmap; it does not alter the
  ## three heterogeneity metrics. Missing CpGs are column-mean imputed solely
  ## for this visualization step.
  clusters <- NULL
  cluster_method <- tolower(as.character(cluster_method[[1L]]))
  k <- suppressWarnings(as.integer(k[[1L]])); if (is.na(k) || k < 2L) k <- 2L
  seed <- suppressWarnings(as.integer(seed[[1L]])); if (is.na(seed)) seed <- 11L
  if (identical(cluster_method, "kmeans") && nrow(meth_mat) >= k) {
    cluster_mat <- meth_mat
    for (jj in seq_len(ncol(cluster_mat))) {
      vv <- cluster_mat[, jj]; mu <- mean(vv, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0.5
      vv[is.na(vv)] <- mu; cluster_mat[, jj] <- vv
    }
    if (nrow(unique(cluster_mat)) >= k) {
      had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      km <- stats::kmeans(cluster_mat, centers = k, nstart = 10L)
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
      clusters <- stats::setNames(as.integer(km$cluster), rownames(meth_mat))
    }
  }
  
  ## Optional per-CpG PDR/qFDRP tables for QC and auditability.
  pdr_by_cpg <- lapply(seq_along(cpg_sites), function(cc) {
    keep <- !is.na(meth_mat[, cc]) & eligible
    vals <- meth_mat[keep, , drop = FALSE]
    disc <- if (nrow(vals)) apply(vals, 1, function(x) {
      x <- x[!is.na(x)]; length(x) >= 4L && length(unique(x)) > 1L
    }) else logical()
    data.frame(Position = cpg_sites[cc], Coverage = sum(weights[keep]),
               PDR = if (any(keep)) sum(weights[keep] * disc) / sum(weights[keep]) else NA_real_)
  }) %>% bind_rows()
  list(
    scores = data.frame(Metric = c("Amplicon PDR", "Window Epipolymorphism", "Amplicon qFDRP"),
                        Value = c(pdr_score, epipoly_score, qfdrp_score)),
    meth_mat = meth_mat, clusters = clusters,
    pdr_by_cpg = pdr_by_cpg, epipolymorphism_by_window = epi_df,
    qfdrp_by_cpg = data.frame()
  )
}

calculate_quma_stats <- function(res_obj, mode = "Sanger") {
  df_long <- res_obj$long_data; df_summary <- res_obj$read_summary %>% filter(!str_detect(Pattern, "^excluded"))
  cpg_sites <- res_obj$genome_info$cpg_pos
  if(nrow(df_summary) == 0) return(list(overall=0, sd_cpg=0, se_cpg=0, sd_seq=0, se_seq=0, cpg_table=data.frame()))
  
  if (mode == "NGS") {
    parsed_count <- suppressWarnings(as.integer(str_extract(
      df_long$ReadID, "(?<=Count)\\d+"
    )))
    parsed_count[is.na(parsed_count) | parsed_count < 1L] <- 1L
    df_long <- df_long %>% mutate(Count = parsed_count)
    overall_meth <- sum(df_long$Methylation * df_long$Count, na.rm = TRUE) /
      sum(df_long$Count, na.rm = TRUE) * 100
    cpg_stats <- df_long %>% group_by(Position) %>% summarise(
      Num_Meth = sum(Methylation * Count, na.rm = TRUE),
      Num_Total = sum(Count, na.rm = TRUE),
      Meth_Pct = Num_Meth / Num_Total * 100,
      .groups = "drop"
    )
  } else {
    overall_meth <- mean(df_long$Methylation) * 100
    cpg_stats <- df_long %>% group_by(Position) %>% summarise(Num_Meth=sum(Methylation), Num_Total=n(), Meth_Pct=mean(Methylation)*100)
  }
  
  sd_cpg <- sd(cpg_stats$Meth_Pct, na.rm=T); se_cpg <- sd_cpg/sqrt(nrow(cpg_stats))
  sd_seq <- sd(df_summary$Meth_Pct, na.rm=T); se_seq <- sd_seq/sqrt(nrow(df_summary))
  full_cpg_stats <- data.frame(Position = cpg_sites) %>% left_join(cpg_stats, by = "Position") %>% replace_na(list(Num_Meth = 0, Num_Total = 0, Meth_Pct = NA))
  return(list(overall=overall_meth, sd_cpg=sd_cpg, se_cpg=se_cpg, sd_seq=sd_seq, se_seq=se_seq, cpg_table=full_cpg_stats))
}

analyze_group_comparison <- function(res_list_A, res_list_B, genome_seq, g1_name="Group 1", g2_name="Group 2") {
  if (is(genome_seq, "DNAStringSet")) genome_seq <- genome_seq[[1]]
  cpg_hits <- matchPattern("CG", genome_seq); cpg_sites <- start(cpg_hits)
  
  agg_long <- function(res_list, group_name) {
    do.call(bind_rows, lapply(names(res_list), function(nm) {
      df <- res_list[[nm]]$long_data
      if(is.null(df) || nrow(df) == 0) return(NULL)
      if(!"Count" %in% names(df)) { cnt <- str_extract(df$ReadID, "(?<=Count)\\d+"); df$Count <- ifelse(is.na(cnt), 1, as.integer(cnt)) }
      df$SampleID <- nm; df$Group <- group_name
      return(df)
    }))
  }
  
  df_A <- agg_long(res_list_A, g1_name); df_B <- agg_long(res_list_B, g2_name)
  if(is.null(df_A) || is.null(df_B)) return(NULL)
  combined_long <- bind_rows(df_A, df_B)
  
  agg_site <- function(df) { df %>% group_by(Position) %>% summarise(Meth = sum(Methylation * Count), Total = sum(Count), Pct = Meth/Total*100) }
  site_A <- agg_site(df_A); site_B <- agg_site(df_B)
  
  full_site_table <- data.frame(Position = cpg_sites) %>%
    left_join(site_A, by="Position") %>% rename(Meth_1=Meth, Total_1=Total, Pct_1=Pct) %>%
    left_join(site_B, by="Position") %>% rename(Meth_2=Meth, Total_2=Total, Pct_2=Pct) %>%
    replace_na(list(Meth_1=0, Total_1=0, Pct_1=0, Meth_2=0, Total_2=0, Pct_2=0))
  
  ## Pooled molecule-level Fisher tests are retained as a fallback and audit
  ## field. When each group has >=2 samples, the displayed site p-values use
  ## sample-level methylation percentages to avoid treating reads as biological
  ## replicates (pseudo-replication).
  sample_site <- function(res_list, group_name) {
    bind_rows(lapply(names(res_list), function(nm) {
      rr <- if (!is.null(res_list[[nm]]$metrics)) res_list[[nm]]$metrics else res_list[[nm]]
      dd <- rr$long_data
      if (is.null(dd) || !nrow(dd)) return(NULL)
      if (!"Count" %in% names(dd)) dd$Count <- 1L
      dd %>% group_by(Position) %>%
        summarise(Pct = sum(Methylation * Count) / sum(Count) * 100,
                  .groups = "drop") %>% mutate(SampleID = nm, Group = group_name)
    }))
  }
  sample_A <- sample_site(res_list_A, g1_name)
  sample_B <- sample_site(res_list_B, g2_name)
  p_vals <- numeric(nrow(full_site_table))
  pooled_p_vals <- numeric(nrow(full_site_table))
  replicate_p_vals <- rep(NA_real_, nrow(full_site_table))
  for(i in 1:nrow(full_site_table)) {
    m1 <- full_site_table$Meth_1[i]; u1 <- full_site_table$Total_1[i] - m1
    m2 <- full_site_table$Meth_2[i]; u2 <- full_site_table$Total_2[i] - m2
    if(full_site_table$Total_1[i] == 0 && full_site_table$Total_2[i] == 0) {
      pooled_p_vals[i] <- NA_real_
    } else {
      mat <- matrix(c(m1, u1, m2, u2), nrow = 2, byrow = TRUE)
      pooled_p_vals[i] <- fisher.test(mat)$p.value
    }
    if (nrow(sample_A) && nrow(sample_B)) {
      a <- sample_A$Pct[sample_A$Position == cpg_sites[i]]
      b <- sample_B$Pct[sample_B$Position == cpg_sites[i]]
      if (sum(is.finite(a)) >= 2L && sum(is.finite(b)) >= 2L)
        replicate_p_vals[i] <- tryCatch(t.test(a, b)$p.value, error = function(e) NA_real_)
    }
  }
  p_vals <- ifelse(is.finite(replicate_p_vals), replicate_p_vals, pooled_p_vals)
  full_site_table$P_Value <- p_vals
  
  # Benjamini-Hochberg FDR correction
  full_site_table$FDR <- p.adjust(p_vals, method = "BH")
  
  get_read_pcts <- function(res_list) {
    unlist(lapply(res_list, function(r) {
      rr <- if (!is.null(r$metrics)) r$metrics else r
      d <- rr$read_summary %>% filter(!str_detect(Pattern, "^excluded"))
      if(nrow(d)==0) return(numeric(0)); return(d$Meth_Pct) 
    }))
  }
  reads_A <- get_read_pcts(res_list_A); reads_B <- get_read_pcts(res_list_B)
  pooled_u_test_p <- NA
  if(length(reads_A) > 0 && length(reads_B) > 0) {
    reads_A <- na.omit(reads_A); reads_B <- na.omit(reads_B)
    if(length(reads_A) > 0 && length(reads_B) > 0) { 
      u_test <- wilcox.test(reads_A, reads_B, exact=FALSE)
      pooled_u_test_p <- u_test$p.value 
    }
  }
  sample_overall <- function(x) {
    vapply(x, function(r) {
      rr <- if (!is.null(r$metrics)) r$metrics else r; d <- rr$read_summary
      d <- d[!grepl("^excluded", d$Pattern) & is.finite(d$Meth_Pct), ]
      if (!nrow(d)) return(NA_real_)
      if ("Count" %in% names(d)) weighted.mean(d$Meth_Pct, d$Count) else mean(d$Meth_Pct)
    }, numeric(1))
  }
  oa <- sample_overall(res_list_A); ob <- sample_overall(res_list_B)
  oa <- oa[is.finite(oa)]; ob <- ob[is.finite(ob)]
  u_test_p <- if (length(oa) >= 2L && length(ob) >= 2L)
    tryCatch(wilcox.test(oa, ob, exact = FALSE)$p.value, error = function(e) NA_real_)
  else pooled_u_test_p
  
  sum_df <- data.frame(
    Group = c(g1_name, g2_name),
    Mean = c(mean(site_A$Pct, na.rm=T), mean(site_B$Pct, na.rm=T)),
    SD = c(sd(oa, na.rm = TRUE), sd(ob, na.rm = TRUE))
  )
  sum_df$SD[!is.finite(sum_df$SD)] <- 0
  
  return(list(site_table = full_site_table, combined_long = combined_long,
              u_test_p = u_test_p, pooled_u_test_p = pooled_u_test_p,
              replicate_p_used = any(is.finite(replicate_p_vals)),
              sample_summary = data.frame(Group = c(g1_name, g2_name),
                                          N = c(length(oa), length(ob)),
                                          Mean = c(mean(oa), mean(ob))),
              summary = sum_df, g1_name=g1_name, g2_name=g2_name))
}

# ==============================================================================
# Plotting Helpers
# ==============================================================================

create_meth_histogram <- function(long_data) {
  if(is.null(long_data) || nrow(long_data) == 0) return(NULL)
  
  df_summary <- long_data %>%
    group_by(ReadID) %>%
    summarise(
      Mean_Meth = mean(Methylation, na.rm = TRUE) * 100,
      Count = if(any(grepl("Count", ReadID))) as.integer(str_extract(first(ReadID), "(?<=Count)\\d+")) else 1
    )
  
  ggplot(df_summary, aes(x=Mean_Meth, weight=Count)) +
    geom_histogram(binwidth=5, fill="steelblue", color="white", alpha=0.9) +
    scale_x_continuous(breaks = seq(0, 100, 10), limits = c(-5, 105)) +
    theme_minimal() +
    labs(
      title = "Methylation Distribution (Per Read)",
      x = "Methylation % per Read",
      y = "Read Count (Frequency)"
    ) +
    theme(text=element_text(size=14))
}

create_abundance_heatmap <- function(long_data) {
  if(is.null(long_data) || nrow(long_data) == 0) return(NULL)
  
  ## Display overview: group retained epialleles by their mean methylation
  ## range.  All reads remain in the analysis; only the visualisation is
  ## aggregated.
  if (!"Count" %in% names(long_data)) {
    long_data$Count <- ifelse(
      grepl("Count", long_data$ReadID),
      as.integer(stringr::str_extract(long_data$ReadID, "(?<=Count)\\d+")),
      1L
    )
  }
  long_data$Count[is.na(long_data$Count) | long_data$Count < 1] <- 1L
  heat_order <- long_data %>%
    group_by(ReadID) %>%
    summarise(Mean_Meth = mean(Methylation, na.rm = TRUE),
              Count = first(Count), .groups = "drop") %>%
    arrange(Mean_Meth, desc(Count), ReadID) %>%
    mutate(ymax = cumsum(Count), ymin = lag(ymax, default = 0)) %>%
    select(ReadID, Count, Mean_Meth, ymin, ymax)
  heat_data <- long_data %>%
    select(ReadID, Position, Methylation) %>%
    left_join(heat_order, by = "ReadID")
  heat_data$Position <- as.numeric(as.character(heat_data$Position))
  cpg_positions <- sort(unique(heat_data$Position[is.finite(heat_data$Position)]))
  x_width <- if (length(cpg_positions) > 1L) min(5, 0.8 * min(diff(cpg_positions))) else 5
  
  p_heat <- ggplot(heat_data) +
    geom_rect(aes(xmin = Position - x_width / 2, xmax = Position + x_width / 2,
                  ymin = ymin, ymax = ymax, fill = factor(Methylation)), color = NA) +
    scale_fill_manual(values = c("0" = "lightblue", "1" = "firebrick"),
                      breaks = c("1", "0"), labels = c("Methylated", "Unmethylated"),
                      name = "Status") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal(base_size = 18) +
    theme(panel.background = element_rect(fill = "grey92", colour = NA),
          panel.grid.major = element_line(colour = "white", linewidth = 0.35),
          panel.grid.minor = element_blank(), axis.text = element_text(size = 14),
          axis.title = element_text(size = 16), axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(
      title = "Abundance heatmap (all retained variants)",
      subtitle = "Band thickness is proportional to read Count; rows ordered by mean methylation",
      x = "CpG Position (bp)", y = "Cumulative read count"
    ) +
    scale_x_continuous(breaks = seq_along(cpg_positions), labels = cpg_positions,
                       expand = c(0, 0))
  p_heat
}

create_lollipop_plot <- function(long_data, point_size = 4, text_size = 14, top_n = 30L) {
  if(is.null(long_data) || nrow(long_data) == 0) return(NULL)
  
  df <- long_data %>%
    mutate(CpG_Index = as.numeric(factor(Position)))
  
  read_stats <- df %>%
    group_by(ReadID) %>%
    summarise(
      Mean_Meth = mean(Methylation, na.rm = TRUE),
      Count = if(any(grepl("Count", ReadID))) as.integer(str_extract(first(ReadID), "(?<=Count)\\d+")) else 1
    ) %>%
    arrange(desc(Count), desc(Mean_Meth), ReadID) %>%
    mutate(Count_Rank = row_number())
  
  top_n <- max(1L, as.integer(top_n))
  ## Select objectively by abundance; methylation order is applied only to
  ## the selected rows for readability.
  keep_ids <- read_stats %>%
    arrange(desc(Count), desc(Mean_Meth), ReadID) %>%
    slice_head(n = top_n) %>% pull(ReadID) %>% as.character()
  df <- df %>% filter(ReadID %in% keep_ids)
  read_stats <- read_stats %>% filter(ReadID %in% keep_ids) %>%
    arrange(Mean_Meth, desc(Count), ReadID)
  
  unique_ids <- read_stats$ReadID
  df$ReadID <- factor(df$ReadID, levels = unique_ids)
  read_stats$ReadID <- factor(read_stats$ReadID, levels = unique_ids)
  
  x_breaks <- unique(df$CpG_Index)
  x_labels <- unique(df$Position)
  
  is_sanger <- all(read_stats$Count == 1)
  
  p1 <- ggplot(df, aes(x=CpG_Index, y=ReadID)) + 
    geom_line(aes(group=ReadID), color="gray80") + 
    geom_point(aes(fill=factor(Methylation)), shape=21, size=point_size, color="black") + 
    scale_fill_manual(values=c("0"="white", "1"="black")) + 
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    theme_minimal() + 
    theme(
      text = element_text(size=text_size), 
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12)
    ) +
    labs(x = "CpG Position (Sequential)", y = if(is_sanger) "Clones / Reads" else "Alleles (methylation state order)")
  
  if (!is_sanger) {
    p1 <- p1 + theme(axis.text.y = element_blank())
    p2 <- ggplot(read_stats, aes(x = Count, y = ReadID)) +
      geom_col(fill = "steelblue") +
      geom_text(aes(label = paste0("#", Count_Rank, " (n=", Count, ")")),
                hjust = -0.1, size = 3.9) +
      theme_void(base_size = 13) + 
      theme(plot.margin = margin(l = 10, r = 45, t = 10, b = 10)) +
      coord_cartesian(clip = "off") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.65))) +
      labs(title = "Read Count")
    return(p1 + p2 + plot_layout(widths = c(3, 1.15)))
  } else {
    return(p1)
  }
}

create_asm_heatmap <- function(het_res) {
  if(is.null(het_res) || is.null(het_res$clusters)) return(NULL)
  mat <- het_res$meth_mat
  if(nrow(mat) < 2) return(NULL)
  cl <- het_res$clusters
  
  row_order <- order(cl, rowMeans(mat, na.rm=TRUE))
  sorted_ids <- rownames(mat)[row_order]
  
  df_plot <- mat %>% 
    rownames_to_column("ReadID") %>% 
    pivot_longer(-ReadID, names_to="Pos", values_to="Meth")
  
  df_plot$Pos <- factor(as.numeric(df_plot$Pos), levels = sort(unique(as.numeric(df_plot$Pos))))
  
  df_plot$Cluster <- cl[df_plot$ReadID]
  df_plot$ReadID <- factor(df_plot$ReadID, levels = sorted_ids)
  
  ggplot(df_plot, aes(x=Pos, y=ReadID, fill=factor(Meth))) + 
    geom_tile(color="white") + 
    scale_fill_manual(values=c("0"="lightblue", "1"="firebrick"), na.value="gray90") + 
    facet_grid(Cluster ~ ., scales="free_y", space="free_y", switch="y") + 
    theme_minimal() + 
    theme(axis.text.y=element_blank(), panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
    labs(fill="Meth", y="Reads (Clustered)", x="CpG Position (Sequential)")
}

create_comp_barplot <- function(summary_df) {
  if(is.null(summary_df) || nrow(summary_df)==0) return(NULL)
  if (!"SD" %in% names(summary_df)) summary_df$SD <- 0
  group_levels <- unique(as.character(summary_df$Group))
  palette_values <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  group_cols <- setNames(rep(palette_values, length.out = length(group_levels)), group_levels)
  ggplot(summary_df, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(width = 0.62, show.legend = FALSE) +
    geom_errorbar(aes(ymin = pmax(0, Mean - SD), ymax = pmin(100, Mean + SD)),
                  width = 0.14, linewidth = 0.5) +
    geom_text(aes(label = formatC(Mean, format = "fg", digits = 3)),
              vjust = -0.45, size = 4) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 100, 20),
                       expand = expansion(mult = c(0, 0.08))) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.margin = margin(12, 12, 18, 12)) +
    labs(y = "Mean methylation (%)", x = NULL,
         title = "Global methylation comparison")
}

create_diff_plot <- function(site_table, g1_name="Group 1", g2_name="Group 2") {
  if(is.null(site_table)) return(NULL)
  df <- site_table %>% mutate(Delta = Pct_1 - Pct_2)
  ggplot(df, aes(x=factor(Position), y=Delta, fill=as.character(Delta>0))) + 
    geom_bar(stat="identity") + 
    scale_fill_manual(
      name = "",
      values = c("TRUE"="firebrick", "FALSE"="steelblue"), 
      breaks = c("TRUE", "FALSE"),
      labels = c(paste(g1_name, "Higher"), paste(g2_name, "Higher"))
    ) + 
    theme_minimal() + 
    labs(title=paste0("Methylation Difference (", g1_name, " vs ", g2_name, ")"), y = "Delta %", x = "CpG Position") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))
}

# ==============================================================================
# 2. UI (Frontend)
# ==============================================================================
ui <- page_navbar(
  title = "PANDA (v1.0.0)", 
  theme = bs_theme(bootswatch = "minty"),
  header = shinyjs::useShinyjs(),
  
  nav_panel(title = "Introduction",
            card(
              h3("PANDA: Phased ANalysis of DNA Amplicons"),
              p(
                em("A unified platform for read-level phased bisulfite analysis of targeted amplicons (Sanger & NGS)"),
                br(), 
                "PANDA preserves methylation patterns at the individual clone or molecule level, rather than reducing an amplicon to only a per-CpG average. The source code, demo data, and detailed documentation are publicly available on ",
                a("https://github.com/kubo-azu/PANDA", 
                  href = "https://github.com/kubo-azu/PANDA", 
                  target = "_blank"),
                "."
              ),
              p(
                strong("Citation: "),
                "Kubota, A., Kobayashi, H., & Tajima, A. (2026). PANDA: Read-Level Phased Analysis of DNA Amplicons for Methylation Studies. ",
                em("bioRxiv"), ". ",
                a("https://doi.org/10.64898/2026.04.01.715790",
                  href = "https://doi.org/10.64898/2026.04.01.715790",
                  target = "_blank")
              ),

              div(class = "row g-3 mb-3",
                  div(class = "col-md-4",
                      div(class = "card h-100 border-primary",
                          div(class = "card-body",
                              h5(class = "card-title", icon("dna"), "Sanger"),
                              p(class = "card-text", "Treats each FASTA record or AB1 base-called sequence as an individual clone. Suitable for phased clonal methylation patterns and small-to-moderate clone sets.")
                          )
                      )
                  ),
                  div(class = "col-md-4",
                      div(class = "card h-100 border-success",
                          div(class = "card-body",
                              h5(class = "card-title", icon("server"), "Amplicon-NGS"),
                              p(class = "card-text", "Dereplicates identical reads, retains their original read counts, and uses those counts for abundance-weighted summaries and plots. The default is to use all retained reads.")
                          )
                      )
                  ),
                  div(class = "col-md-4",
                      div(class = "card h-100 border-info",
                          div(class = "card-body",
                              h5(class = "card-title", icon("chart-line"), "Read-level heterogeneity"),
                              p(class = "card-text", "Reports Amplicon PDR, Window Epipolymorphism, and Amplicon qFDRP to summarize discordant and diverse methylation states across individual molecules.")
                          )
                      )
                  )
              ),
              
              div(class="alert alert-success",
                  h4(icon("shield-halved"), " Data Handling & Privacy"),
                  p(style="font-size: 0.95em;",
                    "Uploaded files are processed in a session-scoped temporary workspace. The application removes uploaded and intermediate files when the session is cleared or terminated, and it does not intentionally retain them as an analysis database. ",
                    strong("This statement describes the application behavior, not the policies of the service hosting it."),
                    " For human or otherwise sensitive data, use the local GUI or CLI under your institution's approved storage and retention policy.")
              ),
              # ----------------------------------------------------------------------
              
              div(class="alert alert-info",
                  h4(icon("flask"), " Quick Start with Demo Data"),
                  p("Download our comprehensive demo dataset (.zip) containing Sanger and NGS samples."),
                  downloadButton("dl_demo_data", "Download Demo Data (PANDA_Demo.zip)", class="btn-primary w-100")
              ),
              
              h4(icon("layer-group"), " Overview"),
              tags$ul(
                tags$li(strong("Sanger Mode:"), " Analysis for clonal sequencing (e.g. TA cloning). Each sequence record is an equally weighted clone (Count = 1)."),
                tags$li(strong("NGS Mode:"), " Optimized for high-depth amplicon sequencing. Identical sequences are dereplicated, while their original read counts are retained as abundance weights."),
                tags$li(strong("Phased outputs:"), " Lollipop plots show methylation patterns in a visibility-oriented CpG order; abundance heatmaps retain physical CpG coordinates and read-count structure."),
                tags$li(strong("Multi-Target & Incremental Upload:"), " Reference supports Multi-FASTA. Files can be accumulated in the queue without resetting."),
                tags$li(strong("Before upload:"), " Trim adapters and primers and perform experiment-specific read QC before analysis. PANDA does not replace a general FASTQ quality-control workflow.")
              ),
              
              h4(icon("list-check"), " How to Use (Step-by-Step)"),
              accordion(
                open = "Step 1: Data Preparation",
                
                accordion_panel("Step 1: Data Preparation",
                                p("The quality of analysis depends on a correct reference amplicon sequence."),
                                div(class="alert alert-warning", icon("triangle-exclamation"), strong(" CRITICAL:"), " Do NOT upload whole genome files (e.g., hg38). Use only amplicon sequences."),
                                
                                h6(strong("Method A: Using Primers information (Recommended)")),
                                tags$ul(
                                  tags$li("Use tools like ", a("UCSC In-Silico PCR", href="https://genome.ucsc.edu/cgi-bin/hgPcr", target="_blank"), "."),
                                  tags$li("Paste your Forward and Reverse primers."),
                                  tags$li("The tool will extract the exact amplicon sequence. Save this as .fasta."),
                                  tags$li("This ensures correct flanking regions and alignment edges.")
                                ),
                                
                                h6(strong("Method B: Manual Extraction")),
                                tags$ol(
                                  tags$li("Search for your target gene (e.g., 'Oct4') in UCSC/NCBI and zoom in to your PCR amplicon region."),
                                  tags$li("Include ", strong("30-50bp flanking regions"), " around your primers to ensure correct alignment edges."),
                                  tags$li("Save as FASTA (.fasta/.fa). Multi-FASTA is supported when multiple reference records are needed.")
                                ),
                                
                                hr(),
                                h6(icon("wand-magic-sparkles"), strong(" Reference Auto-Correction (Strand Detection)")),
                                p("PANDA automatically detects if your sequenced reads are derived from the opposite strand and internally generates the Reverse Complement reference for alignment. (A warning notification will appear if auto-correction is triggered).")
                ),
                
                accordion_panel("Step 2: Prepare Input Reads (Best Practices)",
                                navset_card_underline(
                                  nav_panel(title = "Sanger Sequencing", icon = icon("dna"),
                                            tags$ul(
                                              tags$li(strong("FASTA:"), " .fasta / .fa files containing base-called sequences. (recommended)"),
                                              tags$li(strong(".ab1 files:"), " Upload raw chromatogram files directly.")
                                            )
                                  ),
                                  nav_panel(title = "NGS (Next-Gen)", icon = icon("server"),
                                            p(strong("Accepted file formats: "), " FASTQ (.fastq/.fq/.gz) and FASTA (.fasta/.fa)"),
                                            div(class="alert alert-success",
                                                h6(strong(" Supported Read Modes & Best Practices for Paired-End Data")),
                                                p("To get the cleanest and most accurate results, please choose the correct mode based on your amplicon length:"),
                                                tags$ul(
                                                  tags$li(strong("1. Merged Reads (Highly Recommended):"), " Please merge your Paired-End reads (R1/R2) using external tools (e.g., fastp, FLASH) before uploading. Select 'Single / Merged' mode. This ensures the highest accuracy across the entire amplicon."),
                                                  tags$li(strong("2. Single-End (Read 1 Only) (Recommended):"), " If your amplicon is short and completely covered by Read 1 (e.g., 150bp read for a 170bp amplicon), or if merging is not feasible, ", strong("uploading ONLY the Read 1 file"), " in 'Single / Merged' mode is the best practice. It prevents double-counting and minimizes noise from the lower-quality R2."),
                                                  tags$li(strong("3. Unmerged Pairs (R1 & R2):"), " ", em("Use ONLY for long amplicons with an unsequenced gap in the middle."), " If you use this mode for short, overlapping amplicons, PANDA will double-count the identical molecules, and the higher error rate of R2 might cause the reads to fail the quality filter (resulting in 'No valid results').")
                                                )
                                            )
                                  )
                                )
                ),
                
                accordion_panel("Step 3: Upload & Analyze",
                                navset_card_underline(
                                  nav_panel("Sanger (Single)",
                                            tags$ol(
                                              tags$li("Upload Reference FASTA."),
                                              tags$li("Select 'Single Sample' in the sidebar."),
                                              tags$li("Upload reads (.fasta/.ab1)."),
                                              tags$li("Click ", strong("Run Single Analysis"), "."),
                                              tags$li("Results appear in the tabs.")
                                            )
                                  ),
                                  nav_panel("Sanger (Multi & Comparison)",
                                            h6(icon("wrench"), " Batch Analysis"),
                                            tags$ol(
                                              tags$li("Upload Reference FASTA."),
                                              tags$li("Select 'Multi-Sample Comparison'."),
                                              tags$li("Upload all FASTA files (Control, Test)."),
                                              tags$li("Click ", strong("Run Batch Analysis"), ".")
                                            ),
                                            h6(icon("scale-balanced"), " Group Comparison"),
                                            tags$ol(
                                              tags$li("Go to 'Compare Groups' in sidebar."),
                                              tags$li("Select Group 1 & Group 2 files."),
                                              tags$li("Click ", strong("Run Comparison"), ".")
                                            )
                                  ),
                                  nav_panel("NGS Mode",
                                            h6(icon("wrench"), " Batch Analysis"),
                                            tags$ol(
                                              tags$li("Upload Reference FASTA."),
                                              tags$li("Select your Read Mode ('Merged' or 'Unmerged Pairs')."),
                                              tags$li("Upload your FASTQ/FASTA files to the batch queue."),
                                              tags$li("Click ", strong("Run Batch Analysis"), ".")
                                            ),
                                            h6(icon("scale-balanced"), " Group Comparison"),
                                            tags$ol(
                                              tags$li("Select Group 1 & Group 2 in sidebar."),
                                              tags$li("Click ", strong("Run Comparison"), ".")
                                            )
                                  )
                                )
                ),
                
                accordion_panel("Advanced: Motif Filter (Haplotype Phasing)",
                                p("PANDA includes a powerful filtering tool to perform ", strong("in-silico genotyping"), " and separate alleles based on sequence polymorphisms (SNPs)."),
                                tags$ul(
                                  tags$li(strong("Why use it? "), "If your amplicon contains a known SNP (e.g., Paternal = A, Maternal = G), you can extract reads/clones originating from only one specific allele. This is highly effective for cleanly visualizing Allele-Specific Methylation (ASM) and Loss of Imprinting (LOI)."),
                                  tags$li(strong("How to use: "), "Enter a specific sequence motif in the sidebar. To define a complex haplotype, you can enter ", strong("multiple motifs separated by commas"), " (e.g., 'ATGC, TCGG'). PANDA will only analyze reads that contain ", em("ALL"), " of the specified motifs."),
                                  tags$li(strong("File Upload: "), "For complex haplotype definitions, you can upload a plain text (.txt) or .csv file containing one motif per line.")
                                ),
                                div(class="alert alert-warning", icon("triangle-exclamation"),
                                    strong(" CRITICAL NOTES FOR FILTERING:"),
                                    tags$ul(
                                      tags$li("You must enter the ", strong("Bisulfite-Converted"), " sequence (i.e., unmethylated 'C's should be written as 'T's)."),
                                      tags$li("Using ", strong("C/T SNPs"), " for filtering is ", strong("NOT recommended"), " as they can be confounded with the methylation status itself. It is strongly recommended to use A/G or other non-C/T polymorphisms to ensure accurate phasing.")
                                    )
                                )
                ),
                
                accordion_panel("Step 4: Results & Interpretation",
                                navset_card_underline(
                                  nav_panel("Visualization & Data",
                                            
                                            h5(icon("file"), "1. Sanger (Single & Multi)"),
                                            h6("Detailed analysis for clonal sequences."),
                                            tags$ul(
                                              tags$li(strong("Lollipop Plot:"), " Visualizes methylation status. CpG sites are displayed sequentially (equally spaced) for clear pattern recognition. (Y-axis shows clone names)."),
                                              tags$li(strong("Heterogeneity & ASM:"), " Contains the ASM heatmap and amplicon-level heterogeneity scores."),
                                              tags$li(strong("Sequence Info & Alignments:"), " Quality control metrics and pairwise text alignments to check mismatches directly on the sequence.")
                                            ),
                                            hr(),
                                            
                                            h5(icon("server"), "2. NGS Mode"),
                                            h6("High-depth amplicon analysis optimized for ASM/Imprinting detection."),
                                            tags$ul(
                                              tags$li(strong("Top Sequences (Lollipop):"), " Shows the most abundant unique alleles. CpG sites are ", strong("equally spaced (QUMA-style)"), " for high visibility, combined with a bar chart of read counts. Ideal for presentations."),
                                              tags$li(strong("ASM Profile (Main):"),
                                                      tags$ul(
                                                        tags$li(strong("Heterogeneity & Distribution:"), " Amplicon PDR, window Epipolymorphism, amplicon qFDRP, and a histogram to check bimodality."),
                                                        tags$li(strong("Abundance Heatmap (Weighted):"), " Visualizes allele proportions. Unlike the Lollipop plot, the X-axis reflects the ", strong("actual genomic distance (bp)"), ", allowing you to assess physical read coverage and drop-outs.")
                                                      )
                                              ),
                                              tags$li(strong("Statistics & Alignments:"), " Detailed tables for CpG rates. The Alignments tab shows raw read sequences, perfect for confirming Motif Filter targets.")
                                            ),
                                            hr(),
                                            
                                            h5(icon("scale-balanced"), "3. Group Comparison"),
                                            p("Statistical comparison (Difference Plot, sample-level comparison, and pooled molecule-level Fisher's Exact as an audit fallback) between two groups across Sanger or NGS batches. When at least two samples are available per group, sample-level tests are used to avoid treating sequencing molecules as biological replicates. ", 
                                              strong("P-values for single-CpG comparisons are adjusted for multiple testing using the Benjamini-Hochberg (FDR) method."))
                                  ),
                                  
                                  nav_panel("Heterogeneity Metrics",
                                            p("PANDA reports abundance-weighted, amplicon-adapted heterogeneity metrics based on Scherer et al. (Nucleic Acids Research, 2020). Sanger clones have Count = 1; dereplicated NGS molecules retain their read counts."),
                                            hr(),
                                            
                                            h5(strong("Amplicon PDR (Proportion of Discordant Reads):")),
                                            p("The abundance-weighted percentage of molecules containing at least four observed CpGs and both methylated and unmethylated states. It summarizes within-molecule discordance across the target amplicon."),
                                            hr(),
                                            
                                            h5(strong("Window Epipolymorphism:")),
                                            p("Calculated independently for each consecutive four-CpG window using abundance-weighted epiallele frequencies, then summarized across eligible windows. High values indicate diverse local methylation patterns."),
                                            hr(),
                                            
                                            h5(strong("Amplicon qFDRP:")),
                                            p("Measures the abundance-weighted normalized Hamming distance between molecule pairs over CpGs observed in both molecules. It is an amplicon-level adaptation for long-range phased patterns; a high value indicates diverse molecule-level methylation states.")
                                  )
                                )
                ),
                
                accordion_panel("Technical Note: Alignment & Processing",
                                
                                h6("1. Alignment & Quality Control"),
                                p("PANDA processes data completely within R, eliminating the need for external aligners."),
                                tags$ul(
                                  tags$li(strong("Alignment Strategy:"), " Employs ", strong("Global Pairwise Alignment"), " via ", code("pwalign"), " on ", em("in silico"), " C-to-T converted sequences. ",
                                          "Highly stringent scoring parameters ", strong("(match = 1, mismatch = -3, gap opening = -10, gap extension = -4)"), 
                                          " are applied to prevent artefactual gapping and ensure exact positional coordinate mapping of CpG sites."),
                                  tags$li(strong("Biological Quality Filtering:"), " To seamlessly integrate both NGS and Sanger data, PANDA ", strong("does not use Phred base quality scores"), " for filtering. Instead, reads are strictly filtered by:",
                                          tags$ul(
                                            tags$li("Sequence Identity (Default: \u2265 90%)"),
                                            tags$li("Bisulfite Conversion Rate of non-CpG cytosines (Default: \u2265 95%)")
                                          )
                                  )
                                ),
                                
                                hr(),
                                
                                h6("2. NGS Data Dereplication"),
                                p("In NGS Amplicon Mode, PANDA groups identical sequences to efficiently handle large datasets and identify unique epialleles."),
                                tags$ul(
                                  tags$li(strong("Dereplication"), ": Reads with 100% identical sequences from end to end are collapsed together."),
                                  tags$li(strong("Rank & Count"), ": Ordered by abundance. 'Count' represents the number of raw reads supporting that specific exact sequence."),
                                  tags$li(strong("Sensitivity to Read Boundaries"), ": Because grouping requires a 100% exact match, any variations in read lengths (e.g., untrimmed adapters) will cause identical biological molecules to be split into separate ranks. ", em("(See 'Step 2' for pre-processing best practices to ensure uniform boundaries)."))
                                )
                )
              )
            )
  ),
  
  nav_panel(title = "Sanger Mode",
            layout_sidebar(
              sidebar = sidebar(
                width = 350,
                div(id = "sanger_sidebar_inputs",
                    h5("Input (Sanger)"),
                    fileInput("sanger_genome", "Reference FASTA"),
                    uiOutput("sanger_target_selector_ui"),
                    radioButtons("sanger_mode", "Analysis Type:", choices = c("Single Sample" = "single", "Multi-Sample Comparison" = "multi")),
                    
                    conditionalPanel("input.sanger_mode == 'single'",
                                     radioButtons("sanger_single_fmt", "Format:", c("FASTA"="fasta", ".ab1"="ab1")),
                                     conditionalPanel("input.sanger_single_fmt == 'fasta'", fileInput("sanger_single_fasta", "Reads (FASTA)", multiple=TRUE)),
                                     conditionalPanel("input.sanger_single_fmt == 'ab1'", fileInput("sanger_single_ab1", "Reads (.ab1)", multiple=TRUE),
                                                      sliderInput("ab1_trim_start", "Trim 5'", 0, 100, 20), sliderInput("ab1_trim_end", "Trim 3'", 0, 100, 20)),
                                     actionButton("run_sanger_single", "Run Single Analysis", class="btn-primary w-100"),
                                     br(), br(),
                                     uiOutput("ui_dl_sanger_single")
                    ),
                    
                    conditionalPanel("input.sanger_mode == 'multi'",
                                     p(style="font-size:0.8em; color:gray", "Accumulate files by uploading. Click 'Clear' to reset list."),
                                     fileInput("sanger_multi_files", "Reads (FASTA)", multiple=TRUE),
                                     actionButton("clear_sanger_files", "Clear Files", class="btn-xs btn-secondary", icon=icon("trash")),
                                     verbatimTextOutput("sanger_file_count_txt"),
                                     hr(),
                                     sliderInput("sanger_conv", "Min Conversion (%)", 0, 100, 95), sliderInput("sanger_ident", "Min Identity (%)", 0, 100, 90),
                                     
                                     hr(),
                                     h5("Motif Filter (Optional)"),
                                     p(style="font-size:0.8em; color:gray", "Upload a list of specific sequences (e.g. SNPs/Haplotypes). Only reads containing ALL patterns will be analyzed."),
                                     fileInput("sanger_motif_file", "Upload Motif List (.txt/.csv)", accept = c("text/plain", "text/csv"), placeholder = "No file selected"),
                                     textInput("sanger_motif_text", "Or enter motifs (comma separated):", placeholder = "e.g., ATGC, TCGG"),
                                     hr(),
                                     
                                     actionButton("run_sanger_multi", "Run Batch Analysis", class="btn-primary w-100"),
                                     br(), br(),
                                     uiOutput("ui_dl_sanger_multi"),
                                     
                                     hr(), h5("Compare Groups"),
                                     uiOutput("sanger_group_ui"), 
                                     layout_columns(col_widths = c(6,6),
                                                    textInput("sanger_grp1_name", "G1 Label:", value="Group 1"),
                                                    textInput("sanger_grp2_name", "G2 Label:", value="Group 2")
                                     ),
                                     actionButton("run_sanger_comp", "Run Comparison", class="btn-success w-100"),
                                     br(), br(),
                                     uiOutput("ui_dl_sanger_comp")
                    )
                ),
                hr(), actionButton("reset_sanger", "Reset All", class="btn-warning w-100", icon=icon("refresh")), br()
              ),
              conditionalPanel("input.sanger_mode == 'single'",
                               navset_card_underline(
                                 nav_panel("Lollipop Plot", plotOutput("sanger_plot", height="550px")),
                                 nav_panel("Heterogeneity & ASM", h5("Heterogeneity Scores"), DTOutput("sanger_het_table"), hr(), h5("ASM Clustering (k=2)"), plotOutput("sanger_asm_plot", height="550px")),
                                 nav_panel("Summary & Stats", verbatimTextOutput("sanger_summary"), DTOutput("sanger_cpg_table")),
                                 nav_panel("Sequence Info", DTOutput("sanger_read_table")),
                                 nav_panel("Alignments", verbatimTextOutput("sanger_pairwise_text"))
                               )
              ),
              conditionalPanel("input.sanger_mode == 'multi'",
                               navset_card_underline(
                                 nav_panel("Batch Summary", plotOutput("sanger_batch_qc", height="400px"), DTOutput("sanger_batch_table")),
                                 nav_panel("Individual Details",
                                           uiOutput("sanger_multi_detail_selector"),
                                           navset_card_tab(
                                             nav_panel("Lollipop", plotOutput("sanger_multi_detail_lollipop", height="550px")),
                                             nav_panel("Heterogeneity & ASM", 
                                                       h5("Heterogeneity"), DTOutput("sanger_multi_detail_het_table"), 
                                                       hr(), h5("ASM Clustering"), plotOutput("sanger_multi_detail_asm_plot", height="550px")),
                                             nav_panel("Summary & Stats", verbatimTextOutput("sanger_multi_detail_summary"), DTOutput("sanger_multi_detail_cpg_table")),
                                             nav_panel("Sequence Info", DTOutput("sanger_multi_detail_read_table")),
                                             nav_panel("Alignments", verbatimTextOutput("sanger_multi_detail_pairwise"))
                                           )
                                 ),
                                 nav_panel("Group Comparison", h5("Differential Methylation"), layout_columns(col_widths=c(6,6), plotOutput("sanger_comp_diff"), plotOutput("sanger_comp_bar")), hr(), h5("Statistics"), verbatimTextOutput("sanger_comp_stat"), DTOutput("sanger_comp_site"))
                               )
              )
            )
  ),
  
  nav_panel(title = "NGS Mode",
            layout_sidebar(
              sidebar = sidebar(
                width = 350,
                div(id = "ngs_sidebar_inputs",
                    h5("Input (Batch)"),
                    fileInput("ngs_genome", "Reference FASTA"),
                    uiOutput("ngs_target_selector_ui"),
                    p(style="font-size:0.8em; color:gray", "Accumulate files by uploading. Click 'Clear' to reset list."),
                    
                    radioButtons("ngs_read_mode", "Read Mode:", choices = c("Single / Merged" = "merged", "Unmerged Paired-end (R1 & R2)" = "unmerged")),
                    
                    conditionalPanel("input.ngs_read_mode == 'merged'",
                                     fileInput("ngs_files_merged", "Reads (FASTQ/FASTA)", multiple = TRUE, accept = c(".fastq", ".fq", ".fasta", ".fa", ".gz"))
                    ),
                    conditionalPanel("input.ngs_read_mode == 'unmerged'",
                                     p(style="font-size:0.8em; color:gray", "Upload ALL R1 and R2 files together. PANDA automatically pairs them based on filenames (e.g., '_R1' & '_R2')."),
                                     fileInput("ngs_files_unmerged", "Upload Unmerged Pairs", multiple = TRUE, accept = c(".fastq", ".fq", ".fasta", ".fa", ".gz"))
                    ),
                    
                    actionButton("clear_ngs_files", "Clear Files", class="btn-xs btn-secondary", icon=icon("trash")),
                    verbatimTextOutput("ngs_file_count_txt"),
                    
                    selectInput("ngs_type", "Type", choices = c("FASTQ (fastq/fq.gz)" = "fastq", "FASTA" = "fasta")),
                    sliderInput("ngs_top_n", "Top N for lollipop plot", 10, 100, 30),
                    numericInput("ngs_min_count", "Minimum read count for metrics", value = 1, min = 1, step = 1),
                    numericInput("ngs_workers", "Alignment workers (parallel)", value = min(.panda_workers, 16L), min = 1, max = 16, step = 1),
                    sliderInput("ngs_conv", "Min Conversion", 0, 100, 95), sliderInput("ngs_ident", "Min Identity", 0, 100, 90),
                    
                    hr(),
                    h5("Motif Filter (Optional)"),
                    p(style="font-size:0.8em; color:gray", "Upload a list of specific sequences (e.g. SNPs/Haplotypes). Only reads containing ALL patterns will be analyzed."),
                    fileInput("ngs_motif_file", "Upload Motif List (.txt/.csv)", accept = c("text/plain", "text/csv"), placeholder = "No file selected"),
                    textInput("ngs_motif_text", "Or enter motifs (comma separated):", placeholder = "e.g., ATGC, TCGG"),
                    
                    actionButton("run_ngs", "Run Batch Analysis", class="btn-danger w-100"),
                    br(), br(),
                    uiOutput("ui_dl_ngs_batch"),
                    
                    hr(), h5("Compare Groups"), uiOutput("ngs_group_ui"), 
                    layout_columns(col_widths = c(6,6),
                                   textInput("ngs_grp1_name", "G1 Label:", value="Group 1"),
                                   textInput("ngs_grp2_name", "G2 Label:", value="Group 2")
                    ),
                    actionButton("run_ngs_comp", "Run Comparison", class="btn-success w-100"),
                    br(), br(),
                    uiOutput("ui_dl_ngs_comp")
                ),
                hr(), actionButton("reset_ngs", "Reset All", class="btn-warning w-100", icon=icon("refresh"))
              ),
              navset_card_underline(
                nav_panel("Batch Summary",
                          p(class = "text-muted", "Click the plot for a larger view."),
                          plotOutput("ngs_batch_qc_plot", height="550px", click = "ngs_batch_qc_click"),
                          DTOutput("ngs_batch_table")),
                
                nav_panel("Single Sample Analysis", 
                          uiOutput("ngs_detail_selector"), 
                          navset_card_tab(
                            nav_panel("Statistics", 
                                      div(style = "overflow-y: auto; max-height: 800px;",
                                          h5("CpG Methylation Rates"), DTOutput("ngs_stats_table"),
                                          br(), hr(), br(),
                                          h5("Read Quality Control"), DTOutput("ngs_read_qc"),
                                          br(), br()
                                      )
                            ),
                            nav_panel("Top Sequences", 
                                      h5("Top Abundant Alleles (Ranked)"),
                                      p(class = "text-muted", "Scroll to inspect all selected alleles; click the plot for a larger view."),
                                      div(style = "max-height:700px; overflow:auto;",
                                          plotOutput("ngs_lollipop_plot", height="900px", click = "ngs_lollipop_click"))
                            ),
                            nav_panel("ASM Profile", 
                                      layout_columns(col_widths=c(12, 12),
                                                     card(card_header("1. Heterogeneity Metrics (PDR, Epipoly, qFDRP)"), DTOutput("ngs_detail_het")),
                                                     card(card_header("2. Methylation Distribution (Histogram)"),
                                                          p(class = "text-muted", "Click the plot for a larger view."),
                                                          plotOutput("ngs_dist_plot", height="450px", click = "ngs_dist_click")),
                                                     card(card_header("3. Abundance Heatmap (Weighted)"),
                                                          p(class = "text-muted", "Scroll to inspect the full plot; click the plot for a larger view."),
                                                          div(style = "max-height:700px; overflow:auto;",
                                                              plotOutput("ngs_abund_heatmap", height="1000px", click = "ngs_abund_heatmap_click")))
                                      )
                            ),
                            nav_panel("Alignments", 
                                      verbatimTextOutput("ngs_alignments_text")
                            )
                          )
                ),
                
                nav_panel("Group Comparison", h5("Differential Methylation"),
                          p(class = "text-muted", "Click either plot for a larger view."),
                          layout_columns(col_widths=c(6,6),
                                         plotOutput("ngs_comp_diff_plot", height="500px", click = "ngs_comp_diff_click"),
                                         plotOutput("ngs_comp_bar_plot", height="500px", click = "ngs_comp_bar_click")),
                          hr(), h5("Statistics"), verbatimTextOutput("ngs_comp_stat_txt"), DTOutput("ngs_comp_site_table"))
              )
            )
  )
)

# ==============================================================================
# 3. SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  ## Full-screen, fit-to-viewport modal styling. This deliberately fits the
  ## complete plot in one view; native CLI dimensions remain available from
  ## the PDF export when physical-size reproduction is required.
  modal_fit_css <- tags$style(HTML(
    ".modal-dialog { max-width: calc(100vw - 20px) !important; width: calc(100vw - 20px) !important; margin: 2vh auto; }
     .modal-content { max-height: 96vh; box-sizing: border-box; }
     .modal-body { overflow: hidden !important; padding: 8px; box-sizing: border-box; }
     .panda-modal-plot { width: 100% !important; max-width: 100% !important; height: 76vh !important; overflow: hidden !important; box-sizing: border-box; }
     .panda-modal-plot .shiny-plot-output { width: 100% !important; max-width: 100% !important; overflow: hidden !important; }
     .panda-modal-plot canvas, .panda-modal-plot img { max-width: 100% !important; height: auto !important; }")
  )
  
  # [SECURITY UPDATE] ----------------------------------------------------------
  # Track the path of temporary files uploaded on the server side and ensure they are deleted
  ## This is session bookkeeping, not application state.  Keeping it as a
  ## regular mutable vector avoids reading a reactiveVal from the
  ## onSessionEnded callback, which has no active reactive context.
  session_temp_files <- character()
  
  # アップロードされたファイルをリストに登録するヘルパー関数
  register_temp_files <- function(paths) {
    valid_paths <- paths[!is.na(paths) & paths != ""]
    if(length(valid_paths) > 0) {
      session_temp_files <<- unique(c(session_temp_files, valid_paths))
    }
  }
  
  # A hook that physically and immediately discards temporary files when the user closes a tab (ends session)
  session$onSessionEnded(function() {
    files_to_delete <- session_temp_files
    if(length(files_to_delete) > 0) {
      unlink(files_to_delete, force = TRUE)
    }
  })
  # ----------------------------------------------------------------------------
  
  output$dl_demo_data <- downloadHandler(
    filename = function() { "PANDA_Demo.zip" },
    content = function(file) {
      if(file.exists("PANDA_Demo.zip")) {
        file.copy("PANDA_Demo.zip", file)
      } else {
        showNotification("Demo data file not found in directory.", type="error", duration = NULL)
      }
    }
  )
  
  stored_sanger <- reactiveVal(data.frame())
  
  stored_ngs <- reactiveVal(data.frame(
    name = character(), 
    datapath = character(), 
    datapath_R2 = character(), 
    orig_name_R1 = character(), 
    orig_name_R2 = character(), 
    stringsAsFactors = FALSE
  ))
  
  observeEvent(input$sanger_multi_files, {
    new_df <- input$sanger_multi_files
    register_temp_files(new_df$datapath) # [SECURITY UPDATE] Register the file path
    old_df <- stored_sanger()
    stored_sanger(bind_rows(old_df, new_df))
  })
  
  # Clear Sanger Files
  observeEvent(input$clear_sanger_files, { 
    # [SECURITY UPDATE] When the Clear button is pressed, the registered files are physically deleted
    unlink(stored_sanger()$datapath, force = TRUE)
    
    shinyjs::reset("sanger_multi_files")
    stored_sanger(data.frame()) 
    sanger_multi_batch(NULL)
    sanger_multi_comp(NULL)
  })
  output$sanger_file_count_txt <- renderText({ n <- nrow(stored_sanger()); if(n==0) "No files selected." else paste(n, "files stored.") })
  
  observeEvent(input$ngs_files_merged, {
    new_df <- input$ngs_files_merged
    if(!is.null(new_df)) {
      register_temp_files(new_df$datapath) # [SECURITY UPDATE] Register the file path
      new_df$datapath_R2 <- NA
      new_df$orig_name_R1 <- new_df$name
      new_df$orig_name_R2 <- NA
      
      old_df <- stored_ngs()
      stored_ngs(bind_rows(old_df, new_df[, c("name", "datapath", "datapath_R2", "orig_name_R1", "orig_name_R2")]))
    }
  })
  
  observeEvent(input$ngs_files_unmerged, {
    new_files <- input$ngs_files_unmerged
    if(is.null(new_files)) return()
    
    register_temp_files(new_files$datapath) # [SECURITY UPDATE] Register the file path
    
    processed_files <- new_files %>%
      mutate(
        ReadDir = case_when(
          grepl("_R1|_1(_|\\.)", name, ignore.case = TRUE) ~ "R1",
          grepl("_R2|_2(_|\\.)", name, ignore.case = TRUE) ~ "R2",
          TRUE ~ "Unknown"
        ),
        BaseName = sub("(_R[12]|_[12])(_[^\\.]+)?\\.(fastq|fq|fasta|fa)(\\.gz)?$", "", name, ignore.case = TRUE)
      )
    
    valid_pairs <- processed_files %>%
      filter(ReadDir != "Unknown") %>%
      group_by(BaseName) %>%
      summarise(
        R1_name = name[ReadDir == "R1"][1],
        R2_name = name[ReadDir == "R2"][1],
        R1_path = datapath[ReadDir == "R1"][1],
        R2_path = datapath[ReadDir == "R2"][1],
        .groups = "drop"
      ) %>%
      filter(!is.na(R1_path) & !is.na(R2_path))
    
    paired_files <- c(valid_pairs$R1_name, valid_pairs$R2_name)
    unpaired_files <- processed_files %>% filter(!(name %in% paired_files))
    
    if (nrow(valid_pairs) > 0) {
      new_rows <- data.frame(
        name = valid_pairs$BaseName,
        datapath = valid_pairs$R1_path,
        datapath_R2 = valid_pairs$R2_path,
        orig_name_R1 = valid_pairs$R1_name,
        orig_name_R2 = valid_pairs$R2_name,
        stringsAsFactors = FALSE
      )
      old_df <- stored_ngs()
      stored_ngs(bind_rows(old_df, new_rows))
      
      msg_paired <- paste(sprintf(" [%s]\n  * R1: %s\n  * R2: %s", valid_pairs$BaseName, valid_pairs$R1_name, valid_pairs$R2_name), collapse="\n\n")
      
      msg_unpaired <- ""
      if (nrow(unpaired_files) > 0) {
        msg_unpaired <- paste("\n\n⚠️ [Unpaired / Ignored Files]\n  *", paste(unpaired_files$name, collapse="\n  * "))
      }
      
      showModal(modalDialog(
        title = "Auto-Pairing Successful!",
        pre(paste0("Successfully created ", nrow(valid_pairs), " pair(s):\n\n", msg_paired, msg_unpaired)),
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      
    } else {
      showNotification("No valid R1/R2 pairs found. Please check filenames.", type = "warning")
    }
  })
  
  # Clear NGS Files
  observeEvent(input$clear_ngs_files, { 
    # [SECURITY UPDATE] When the Clear button is pressed, the registered files are physically deleted
    unlink(c(stored_ngs()$datapath, stored_ngs()$datapath_R2), force = TRUE)
    
    shinyjs::reset("ngs_files_merged")
    shinyjs::reset("ngs_files_unmerged")
    stored_ngs(data.frame(
      name = character(), 
      datapath = character(), 
      datapath_R2 = character(), 
      orig_name_R1 = character(), 
      orig_name_R2 = character(), 
      stringsAsFactors = FALSE
    )) 
    ngs_batch_results(NULL)
    ngs_comp_results(NULL)
  })
  output$ngs_file_count_txt <- renderText({ n <- nrow(stored_ngs()); if(n==0) "No files selected." else paste(n, "files stored.") })
  
  sanger_single_res <- reactiveVal(NULL); sanger_single_het <- reactiveVal(NULL)
  sanger_multi_batch <- reactiveVal(NULL); sanger_multi_comp <- reactiveVal(NULL)
  
  sanger_genomes_list <- reactiveVal(NULL)
  observeEvent(input$sanger_genome, {
    if (!is.null(input$sanger_genome)) {
      register_temp_files(input$sanger_genome$datapath) # [SECURITY UPDATE]
      sanger_genomes_list(readDNAStringSet(input$sanger_genome$datapath))
      # [SECURITY UPDATE] Files can be deleted immediately after being loaded into memory
      unlink(input$sanger_genome$datapath, force = TRUE) 
    }
  })
  
  output$sanger_target_selector_ui <- renderUI({ 
    gl <- sanger_genomes_list()
    if(is.null(gl)) return(NULL)
    g_names <- names(gl)
    if(is.null(g_names)) g_names <- paste0("Target_", seq_along(gl))
    selectInput("sanger_target_select", "Select Target Amplicon:", choices = g_names) 
  })
  
  get_sanger_target_seq <- reactive({ 
    req(sanger_genomes_list(), input$sanger_target_select)
    all_g <- sanger_genomes_list()
    if(is.null(names(all_g))) names(all_g) <- paste0("Target_", seq_along(all_g))
    all_g[[input$sanger_target_select]] 
  })
  
  # Reset Sanger
  observeEvent(input$reset_sanger, { 
    # [SECURITY UPDATE] Files are physically deleted during reset
    unlink(stored_sanger()$datapath, force = TRUE)
    
    shinyjs::reset("sanger_sidebar_inputs")
    stored_sanger(data.frame())
    sanger_single_res(NULL)
    sanger_single_het(NULL)
    sanger_multi_batch(NULL)
    sanger_multi_comp(NULL)
    sanger_genomes_list(NULL)
    showNotification("Sanger Mode Reset Complete", type="warning") 
  })
  
  # Run Sanger Single
  observeEvent(input$run_sanger_single, {
    req(get_sanger_target_seq()); withProgress(message = 'Processing...', value = 0, {
      tryCatch({
        genome <- get_sanger_target_seq(); incProgress(0.3, detail = "Reading...")
        if(input$sanger_single_fmt=="fasta") { 
          req(input$sanger_single_fasta); 
          register_temp_files(input$sanger_single_fasta$datapath) # [SECURITY UPDATE]
          reads <- readDNAStringSet(input$sanger_single_fasta$datapath) 
        } else { 
          req(input$sanger_single_ab1); 
          register_temp_files(input$sanger_single_ab1$datapath) # [SECURITY UPDATE]
          reads <- process_ab1_files(input$sanger_single_ab1$datapath, input$sanger_single_ab1$name, input$ab1_trim_start, input$ab1_trim_end) 
        }
        incProgress(0.6, detail = "Aligning..."); res <- run_bisulfite_alignment(genome, reads, input$sanger_ident, input$sanger_conv, return_alignments = TRUE) 
        if(isTRUE(res$was_flipped)) showNotification("Reference automatically flipped to Reverse Complement.", type="warning", duration=10)
        
        sanger_single_res(res); sanger_single_het(calculate_heterogeneity(res)); incProgress(1, detail = "Done!")
        
      }, error = function(e) showNotification(e$message, type="error", duration = NULL))
    })
  })
  
  output$sanger_summary <- renderPrint({ req(sanger_single_res()); stats <- calculate_quma_stats(sanger_single_res(), "Sanger"); cat(sprintf("Overall Methylation: %.1f%%\n", stats$overall)) })
  output$sanger_plot <- renderPlot({ req(sanger_single_res()); create_lollipop_plot(sanger_single_res()$long_data) })
  
  output$sanger_het_table <- renderDT({ req(sanger_single_het()); datatable(sanger_single_het()$scores, options=list(dom='t'), rownames=FALSE) %>% formatRound(2, 3) })
  
  output$sanger_asm_plot <- renderPlot({ req(sanger_single_het()); create_asm_heatmap(sanger_single_het()) })
  output$sanger_cpg_table <- renderDT({ req(sanger_single_res()); datatable(calculate_quma_stats(sanger_single_res(), "Sanger")$cpg_table, rownames=FALSE) %>% formatRound("Meth_Pct", 1) })
  output$sanger_read_table <- renderDT({ req(sanger_single_res()); datatable(sanger_single_res()$read_summary, options=list(scrollX=T), rownames=FALSE) })
  output$sanger_pairwise_text <- renderPrint({ req(sanger_single_res()); cat(paste(sanger_single_res()$alignments$pairwise, collapse="\n----------------\n")) })
  
  # Run Sanger Multi Batch
  observeEvent(input$run_sanger_multi, {
    req(get_sanger_target_seq()); files <- stored_sanger(); if(nrow(files)==0) { showNotification("No files selected.", type="error", duration = NULL); return() }
    genome <- get_sanger_target_seq(); res_list <- list(); het_list <- list(); sum_list <- list()
    
    target_motifs <- character()
    if (!is.null(input$sanger_motif_file)) {
      tryCatch({
        register_temp_files(input$sanger_motif_file$datapath) # [SECURITY UPDATE]
        target_motifs <- c(target_motifs, readLines(input$sanger_motif_file$datapath))
        unlink(input$sanger_motif_file$datapath, force = TRUE) # メモリ展開後に即時削除
      }, error = function(e) showNotification("Error reading motif file", type="warning"))
    }
    if (input$sanger_motif_text != "") {
      target_motifs <- c(target_motifs, strsplit(input$sanger_motif_text, ",")[[1]])
    }
    target_motifs <- trimws(target_motifs)
    target_motifs <- target_motifs[target_motifs != ""]
    withProgress(message = paste('Batch Analysis:', input$sanger_target_select), value = 0, {
      n <- nrow(files); for(i in 1:n) {
        incProgress(1/n, detail = files$name[i])
        try({
          reads <- readDNAStringSet(files$datapath[i])
          
          if (length(target_motifs) > 0) {
            reads_char <- as.character(reads)
            keep_idx <- sapply(reads_char, function(r) {
              all(sapply(target_motifs, function(m) grepl(m, r)))
            })
            reads <- reads[keep_idx]
          }
          
          if (length(reads) > 0) {
            res <- run_bisulfite_alignment(genome, reads, input$sanger_ident, input$sanger_conv, return_alignments = TRUE) 
            if(!is.null(res)) { 
              nm <- files$name[i]; res_list[[nm]] <- res
              het <- calculate_heterogeneity(res); het_list[[nm]] <- het
              stats <- calculate_quma_stats(res, "Sanger")
              
              sum_list[[length(sum_list)+1]] <- data.frame(
                Sample=nm, 
                Total=res$counts$total, 
                Used=res$counts$used, 
                Meth=round(stats$overall, 1),
                PDR=round(het_list[[nm]]$scores$Value[1], 2),
                Epipoly=round(het_list[[nm]]$scores$Value[2], 4),
                qFDRP=round(het_list[[nm]]$scores$Value[3], 4)
              ) 
              if(isTRUE(res$was_flipped)) showNotification(paste(nm, ": Reference automatically flipped to Reverse Complement."), type="warning", duration=10)
            }
          }
        })
      }
    }); 
    
    if(length(res_list)>0) { 
      sanger_multi_batch(list(details=res_list, hets=het_list, summary=bind_rows(sum_list), genome=genome)) 
      if(length(target_motifs) > 0) showNotification(paste("Filtered by", length(target_motifs), "motifs."), type="message")
    } else { 
      showNotification("No valid results. Check file format or filter criteria.", type="error", duration = NULL) 
    }
  })
  
  output$sanger_batch_qc <- renderPlot({
    req(sanger_multi_batch())
    d <- sanger_multi_batch()$summary %>% filter(is.finite(Meth))
    validate(need(nrow(d) > 0, "No finite methylation values are available for plotting."))
    ggplot(d, aes(x=Sample, y=Meth, fill=Sample)) +
      geom_col(colour="black", show.legend=FALSE) +
      scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20)) + theme_minimal(base_size=14)
  })
  output$sanger_batch_table <- renderDT({ req(sanger_multi_batch()); datatable(sanger_multi_batch()$summary, rownames=FALSE) })
  
  output$sanger_multi_detail_selector <- renderUI({ req(sanger_multi_batch()); selectInput("sanger_multi_view_sample", "Select Sample:", choices=names(sanger_multi_batch()$details)) })
  get_sanger_multi_sel <- reactive({ req(sanger_multi_batch(), input$sanger_multi_view_sample); list(res = sanger_multi_batch()$details[[input$sanger_multi_view_sample]], het = sanger_multi_batch()$hets[[input$sanger_multi_view_sample]]) })
  output$sanger_multi_detail_lollipop <- renderPlot({ req(get_sanger_multi_sel()); create_lollipop_plot(get_sanger_multi_sel()$res$long_data) })
  
  output$sanger_multi_detail_het_table <- renderDT({ datatable(get_sanger_multi_sel()$het$scores, options=list(dom='t'), rownames=FALSE) %>% formatRound(2, 3) })
  
  output$sanger_multi_detail_asm_plot <- renderPlot({ req(get_sanger_multi_sel()); create_asm_heatmap(get_sanger_multi_sel()$het) })
  output$sanger_multi_detail_summary <- renderPrint({ stats <- calculate_quma_stats(get_sanger_multi_sel()$res, "Sanger"); cat(sprintf("Overall Methylation: %.1f%%\n", stats$overall)) })
  output$sanger_multi_detail_cpg_table <- renderDT({ datatable(calculate_quma_stats(get_sanger_multi_sel()$res, "Sanger")$cpg_table, rownames=FALSE) %>% formatRound("Meth_Pct", 1) })
  output$sanger_multi_detail_read_table <- renderDT({ datatable(get_sanger_multi_sel()$res$read_summary, options=list(scrollX=T), rownames=FALSE) })
  output$sanger_multi_detail_pairwise <- renderPrint({ cat(paste(get_sanger_multi_sel()$res$alignments$pairwise, collapse="\n----------------\n")) })
  
  output$sanger_group_ui <- renderUI({ req(sanger_multi_batch()); samps <- names(sanger_multi_batch()$details); tagList(selectInput("sanger_grp1", "Group 1:", samps, multiple=T), selectInput("sanger_grp2", "Group 2:", samps, multiple=T)) })
  
  # Run Sanger Comp
  observeEvent(input$run_sanger_comp, { 
    if (is.null(sanger_multi_batch())) {
      showNotification("Error: Please run Batch Analysis first.", type = "error", duration = NULL)
      return()
    }
    if (is.null(input$sanger_grp1) || is.null(input$sanger_grp2)) {
      showNotification("Warning: Please select at least one sample for both Group 1 and Group 2.", type = "warning")
      return()
    }
    
    withProgress(message = 'Comparing...', value = 0.5, { 
      sanger_multi_comp(analyze_group_comparison(
        sanger_multi_batch()$details[input$sanger_grp1], 
        sanger_multi_batch()$details[input$sanger_grp2], 
        sanger_multi_batch()$genome,
        g1_name = input$sanger_grp1_name, 
        g2_name = input$sanger_grp2_name
      )) 
    }) 
  })
  
  output$sanger_comp_diff <- renderPlot({ req(sanger_multi_comp()); create_diff_plot(sanger_multi_comp()$site_table, sanger_multi_comp()$g1_name, sanger_multi_comp()$g2_name) })
  output$sanger_comp_bar <- renderPlot({ req(sanger_multi_comp()); create_comp_barplot(sanger_multi_comp()$summary) })
  output$sanger_comp_stat <- renderPrint({ 
    req(sanger_multi_comp())
    s <- sanger_multi_comp()$summary
    cat(sprintf("%s Mean: %.1f%%, %s Mean: %.1f%%\nMann-Whitney P-val: %g", 
                s$Group[1], s$Mean[1], s$Group[2], s$Mean[2], sanger_multi_comp()$u_test_p)) 
  })
  
  output$sanger_comp_site <- renderDT({ 
    req(sanger_multi_comp())
    comp <- sanger_multi_comp()
    df <- comp$site_table
    colnames(df) <- c("Position", paste(comp$g1_name, "Meth"), paste(comp$g1_name, "Total"), paste(comp$g1_name, "Pct"), paste(comp$g2_name, "Meth"), paste(comp$g2_name, "Total"), paste(comp$g2_name, "Pct"), "P_Value", "FDR")
    datatable(df, rownames=FALSE) %>% formatRound(c(4, 7), 1) %>% formatSignif(c(8, 9), 3) 
  })
  
  ngs_batch_results <- reactiveVal(NULL); ngs_comp_results <- reactiveVal(NULL)
  
  ngs_genomes_list <- reactiveVal(NULL)
  observeEvent(input$ngs_genome, {
    if (!is.null(input$ngs_genome)) {
      register_temp_files(input$ngs_genome$datapath) # [SECURITY UPDATE]
      ngs_genomes_list(readDNAStringSet(input$ngs_genome$datapath))
      unlink(input$ngs_genome$datapath, force = TRUE) # メモリ展開後に即時削除
    }
  })
  
  output$ngs_target_selector_ui <- renderUI({ 
    gl <- ngs_genomes_list()
    if(is.null(gl)) return(NULL)
    g_names <- names(gl)
    if(is.null(g_names)) g_names <- paste0("Target_", seq_along(gl))
    selectInput("ngs_target_select", "Select Target Amplicon:", choices = g_names) 
  })
  
  # Reset NGS
  observeEvent(input$reset_ngs, { 
    # [SECURITY UPDATE] リセット時にもファイルを物理的に削除
    unlink(c(stored_ngs()$datapath, stored_ngs()$datapath_R2), force = TRUE)
    
    shinyjs::reset("ngs_sidebar_inputs")
    stored_ngs(data.frame(
      name = character(), 
      datapath = character(), 
      datapath_R2 = character(), 
      orig_name_R1 = character(), 
      orig_name_R2 = character(), 
      stringsAsFactors = FALSE
    ))
    ngs_batch_results(NULL)
    ngs_comp_results(NULL)
    ngs_genomes_list(NULL)
    showNotification("NGS Mode Reset Complete", type="warning") 
  })
  
  # Run NGS Batch
  observeEvent(input$run_ngs, {
    req(ngs_genomes_list(), input$ngs_target_select); files <- stored_ngs(); if(nrow(files)==0) { showNotification("No files.", type="error", duration = NULL); return() }
    
    all_genomes <- ngs_genomes_list()
    if(is.null(names(all_genomes))) names(all_genomes) <- paste0("Target_", seq_along(all_genomes))
    target_genome <- all_genomes[[input$ngs_target_select]]
    
    target_motifs <- character()
    if (!is.null(input$ngs_motif_file)) {
      tryCatch({
        register_temp_files(input$ngs_motif_file$datapath) # [SECURITY UPDATE]
        target_motifs <- c(target_motifs, readLines(input$ngs_motif_file$datapath))
        unlink(input$ngs_motif_file$datapath, force = TRUE) # メモリ展開後に即時削除
      }, error = function(e) showNotification("Error reading motif file", type="warning"))
    }
    if (input$ngs_motif_text != "") {
      target_motifs <- c(target_motifs, strsplit(input$ngs_motif_text, ",")[[1]])
    }
    target_motifs <- trimws(target_motifs)
    target_motifs <- target_motifs[target_motifs != ""]
    
    ## This value belongs to the NGS observer. Defining it here prevents the
    ## batch loop from referring to a variable created in the Sanger scope.
    min_count <- suppressWarnings(as.integer(input$ngs_min_count))
    if (is.na(min_count) || min_count < 1L) min_count <- 1L
    workers <- suppressWarnings(as.integer(input$ngs_workers))
    if (is.na(workers) || workers < 1L) workers <- .panda_workers
    workers <- min(workers, 16L)
    
    res_list <- list(); het_list <- list(); sum_list <- list()
    
    withProgress(message = 'NGS Batch Analysis...', value = 0, {
      n <- nrow(files)
      for(i in 1:n) {
        incProgress(1/n, detail = files$name[i])
        
        tryCatch({
          reads_char <- .panda_read_sequences(
            files$datapath[i], files$orig_name_R1[i], input$ngs_type
          )
          is_unmerged_paired <- (input$ngs_read_mode == "unmerged" && !is.na(files$datapath_R2[i]))
          all_reads_set <- DNAStringSet()
          
          if (is_unmerged_paired) {
            reads_char_r2 <- .panda_read_sequences(
              files$datapath_R2[i], files$orig_name_R2[i], input$ngs_type
            )
            min_len <- min(length(reads_char), length(reads_char_r2))
            
            if (min_len > 0) {
              reads_char <- reads_char[1:min_len]
              reads_char_r2 <- reads_char_r2[1:min_len]
              
              r2_set <- DNAStringSet(reads_char_r2)
              r2_rc <- as.character(reverseComplement(r2_set))
              
              paired_seqs <- paste(reads_char, r2_rc, sep="---PAIR---")
              
              if (length(target_motifs) > 0) {
                keep_idx <- sapply(paired_seqs, function(r) {
                  all(sapply(target_motifs, function(m) grepl(m, r)))
                })
                paired_seqs <- paired_seqs[keep_idx]
              }
              
              if (length(paired_seqs) > 0) {
                counts <- sort(table(paired_seqs), decreasing = TRUE)
                counts <- counts[counts >= min_count]
                all_reads_set <- .panda_expand_dereplicated(counts, paired = TRUE)
                top_seqs <- head(counts, input$ngs_top_n)
                
                expanded_names <- character()
                expanded_seqs <- character()
                
                for(seq_idx in seq_along(top_seqs)) {
                  pair_str <- names(top_seqs)[seq_idx]
                  cnt <- as.integer(top_seqs[seq_idx])
                  base_name <- paste0("Rank", seq_idx, "_Count", cnt)
                  
                  parts <- strsplit(pair_str, "---PAIR---")[[1]]
                  
                  expanded_names <- c(expanded_names, paste0(base_name, "_R1"), paste0(base_name, "_R2"))
                  expanded_seqs <- c(expanded_seqs, parts[1], parts[2])
                }
                
                reads_set <- DNAStringSet(expanded_seqs)
                names(reads_set) <- expanded_names
              } else {
                all_reads_set <- DNAStringSet()
                reads_set <- DNAStringSet()
              }
            } else {
              all_reads_set <- DNAStringSet()
              reads_set <- DNAStringSet()
            }
          } else {
            if (!is.na(files$datapath_R2[i])) {
              reads_char_r2 <- .panda_read_sequences(
                files$datapath_R2[i], files$orig_name_R2[i], input$ngs_type
              )
              r2_set <- DNAStringSet(reads_char_r2)
              r2_rc <- as.character(reverseComplement(r2_set))
              reads_char <- c(reads_char, r2_rc) 
            }
            
            if (length(target_motifs) > 0) {
              keep_idx <- sapply(reads_char, function(r) {
                all(sapply(target_motifs, function(m) grepl(m, r)))
              })
              reads_char <- reads_char[keep_idx]
            }
            
            if (length(reads_char) > 0) {
              counts <- sort(table(reads_char), decreasing = TRUE)
              counts <- counts[counts >= min_count]
              all_reads_set <- .panda_expand_dereplicated(counts)
              reads_set <- .panda_expand_dereplicated(counts, input$ngs_top_n)
            } else {
              all_reads_set <- DNAStringSet()
              reads_set <- DNAStringSet()
            }
          }
          
          if (length(reads_set) > 0) {
            ## Align the complete dereplicated set once.  The previous code
            ## aligned Top-N and then aligned all sequences again for metrics.
            ## The shared alignment cache prevents repeated work across files,
            ## while this single full pass removes the within-sample duplicate.
            res <- run_bisulfite_alignment(
              target_genome,
              all_reads_set,
              input$ngs_ident,
              input$ngs_conv,
              return_alignments = TRUE,
              workers = workers
            )
            
            if(!is.null(res)) {
              if(isTRUE(res$was_flipped)) showNotification(paste(files$name[i], ": Reference automatically flipped to Reverse Complement."), type="warning", duration=10)
              
              if (is_unmerged_paired) {
                res$read_summary <- res$read_summary %>%
                  mutate(BaseID = sub("_R[12]$", "", ReadID)) %>%
                  group_by(BaseID) %>%
                  summarise(
                    Strand = paste(unique(Strand), collapse="/"),
                    Mismatches = sum(Mismatches),
                    Gaps = sum(Gaps),
                    Identity_Pct = mean(Identity_Pct),
                    Meth_Pct = mean(Meth_Pct, na.rm=TRUE),
                    Conv_Pct = mean(Conv_Pct),
                    CpG_Count = sum(CpG_Count),
                    Pattern = if(any(grepl("excluded", Pattern))) "excluded (Pair Failed QC)" else "Passed",
                    .groups = "drop"
                  ) %>%
                  rename(ReadID = BaseID)
                
                res$long_data <- res$long_data %>%
                  mutate(ReadID = sub("_R[12]$", "", ReadID))
                
                excluded_ids <- res$read_summary$ReadID[grepl("excluded", res$read_summary$Pattern)]
                res$long_data <- res$long_data %>%
                  filter(!ReadID %in% excluded_ids)
                
                res$counts$total <- res$counts$total / 2
              }
              
              res$long_data <- res$long_data %>% mutate(Count = as.integer(str_extract(ReadID, "(?<=Count)\\d+")))
              res$read_summary <- res$read_summary %>% mutate(Count = as.integer(str_extract(ReadID, "(?<=Count)\\d+")))
              
              ## Quantitative summaries use the complete dereplicated result.
              metrics_res <- res
              if (is.null(metrics_res)) metrics_res <- res
              
              ## Keep all retained variants for distribution and heatmap plots.
              ## Top-N is applied only inside the lollipop plot.
              display_res <- res
              
              nm <- files$name[i]
              display_res$metrics <- metrics_res
              res_list[[nm]] <- display_res
              het_list[[nm]] <- calculate_heterogeneity(metrics_res)
              stats <- calculate_quma_stats(metrics_res, "NGS")
              
              real_total_reads <- sum(metrics_res$read_summary$Count)
              real_used_reads <- sum(metrics_res$read_summary$Count[metrics_res$read_summary$Pattern == "Passed" | !grepl("^excluded", metrics_res$read_summary$Pattern)])
              
              sum_list[[length(sum_list)+1]] <- data.frame(
                Sample=nm, 
                Total=real_total_reads, 
                Used=real_used_reads, 
                Meth=round(stats$overall, 1), 
                PDR=round(het_list[[nm]]$scores$Value[1], 2),
                Epipoly=round(het_list[[nm]]$scores$Value[2], 4),
                qFDRP=round(het_list[[nm]]$scores$Value[3], 4)
              )
            }
          }
          gc() 
          
        }, error = function(e) { 
          showNotification(paste("Error in", files$name[i], ":", e$message), type="error", duration = NULL) 
        })
      }
    })
    
    if(length(res_list)>0) { 
      ngs_batch_results(list(details=res_list, hets=het_list, summary=bind_rows(sum_list), genome=target_genome)) 
      if(length(target_motifs) > 0) showNotification(paste("Filtered by", length(target_motifs), "motifs."), type="message")
    } else { 
      showNotification("No valid results. Check file format or filter criteria.", type="error", duration = NULL) 
    }
  })
  
  output$ngs_batch_qc_plot <- renderPlot({
    req(ngs_batch_results())
    summary_plot <- ngs_batch_results()$summary %>% filter(is.finite(Meth))
    validate(need(nrow(summary_plot) > 0, "No finite methylation values are available for plotting."))
    ggplot(summary_plot, aes(x=Sample, y=Meth, fill=Sample)) +
      geom_col(colour="black", show.legend=FALSE) +
      scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20)) +
      theme_minimal(base_size=14) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
  })
  output$ngs_batch_table <- renderDT({ req(ngs_batch_results()); datatable(ngs_batch_results()$summary, rownames=FALSE) })
  output$ngs_detail_selector <- renderUI({ req(ngs_batch_results()); selectInput("ngs_view_sample", "Select Sample:", choices=names(ngs_batch_results()$details)) })
  get_ngs_sel <- reactive({ req(ngs_batch_results(), input$ngs_view_sample); list(res = ngs_batch_results()$details[[input$ngs_view_sample]], het = ngs_batch_results()$hets[[input$ngs_view_sample]]) })
  
  output$ngs_dist_plot <- renderPlot({ req(get_ngs_sel()); create_meth_histogram(get_ngs_sel()$res$long_data) })
  output$ngs_abund_heatmap <- renderPlot({ req(get_ngs_sel()); create_abundance_heatmap(get_ngs_sel()$res$long_data) })
  output$ngs_lollipop_plot <- renderPlot({ req(get_ngs_sel()); create_lollipop_plot(get_ngs_sel()$res$long_data, top_n = input$ngs_top_n) })
  
  ## A click opens a larger modal without changing the underlying data or
  ## ordering. The in-card plot remains scrollable for routine inspection.
  observeEvent(input$ngs_abund_heatmap_click, {
    showModal(modalDialog(
      title = "Abundance heatmap (full size)",
      modal_fit_css,
      div(class = "panda-modal-plot",
          plotOutput("ngs_abund_heatmap_modal", width = "100%", height = "76vh")),
      size = "l", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$ngs_abund_heatmap_modal <- renderPlot({
    req(get_ngs_sel())
    create_abundance_heatmap(get_ngs_sel()$res$long_data)
  }, res = 110)
  
  observeEvent(input$ngs_lollipop_click, {
    showModal(modalDialog(
      title = "Top abundant alleles (full size)",
      modal_fit_css,
      div(class = "panda-modal-plot",
          plotOutput("ngs_lollipop_modal", width = "100%", height = "76vh")),
      size = "l", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$ngs_lollipop_modal <- renderPlot({
    req(get_ngs_sel())
    create_lollipop_plot(get_ngs_sel()$res$long_data, top_n = input$ngs_top_n)
  }, res = 110)
  
  observeEvent(input$ngs_batch_qc_click, {
    showModal(modalDialog(
      title = "NGS batch summary (full size)",
      modal_fit_css,
      div(class = "panda-modal-plot",
          plotOutput("ngs_batch_qc_modal", width = "100%", height = "76vh")),
      size = "l", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$ngs_batch_qc_modal <- renderPlot({
    req(ngs_batch_results())
    d <- ngs_batch_results()$summary %>% filter(is.finite(Meth))
    validate(need(nrow(d) > 0, "No finite methylation values are available for plotting."))
    ggplot(d, aes(x = Sample, y = Meth, fill = Sample)) +
      geom_col(colour = "black", show.legend = FALSE) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      theme_minimal(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13)) +
      labs(x = "Sample", y = "Methylation (%)")
  }, res = 110)
  
  observeEvent(input$ngs_dist_click, {
    showModal(modalDialog(
      title = "Methylation distribution (full size)",
      modal_fit_css,
      div(class = "panda-modal-plot",
          plotOutput("ngs_dist_modal", width = "100%", height = "76vh")),
      size = "l", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$ngs_dist_modal <- renderPlot({
    req(get_ngs_sel())
    create_meth_histogram(get_ngs_sel()$res$long_data)
  }, res = 110)
  
  output$ngs_detail_het <- renderDT({ datatable(get_ngs_sel()$het$scores, options=list(dom='t'), rownames=FALSE) %>% formatRound(2, 3) })
  
  output$ngs_stats_table <- renderDT({ datatable(calculate_quma_stats(get_ngs_sel()$res, "NGS")$cpg_table, rownames=FALSE) %>% formatRound("Meth_Pct", 1) })
  output$ngs_read_qc <- renderDT({ datatable(get_ngs_sel()$res$read_summary, options=list(scrollX=T), rownames=FALSE) })
  
  output$ngs_alignments_text <- renderPrint({ req(get_ngs_sel()); cat(paste(get_ngs_sel()$res$alignments$pairwise, collapse="\n----------------\n")) })
  
  output$ngs_group_ui <- renderUI({ req(ngs_batch_results()); tagList(selectInput("ngs_grp1", "G1:", names(ngs_batch_results()$details), multiple=T), selectInput("ngs_grp2", "G2:", names(ngs_batch_results()$details), multiple=T)) })
  
  # Run NGS Comp
  observeEvent(input$run_ngs_comp, { 
    if (is.null(ngs_batch_results())) {
      showNotification("Error: Please run Batch Analysis first.", type = "error", duration = NULL)
      return()
    }
    if (is.null(input$ngs_grp1) || is.null(input$ngs_grp2)) {
      showNotification("Warning: Please select at least one sample for both Group 1 and Group 2.", type = "warning")
      return()
    }
    
    withProgress(message = 'Comparing...', value = 0.5, {
      ngs_comp_results(analyze_group_comparison(
        ngs_batch_results()$details[input$ngs_grp1], 
        ngs_batch_results()$details[input$ngs_grp2], 
        ngs_batch_results()$genome,
        g1_name = input$ngs_grp1_name,
        g2_name = input$ngs_grp2_name
      )) 
    })
  })
  
  output$ngs_comp_diff_plot <- renderPlot({ req(ngs_comp_results()); create_diff_plot(ngs_comp_results()$site_table, ngs_comp_results()$g1_name, ngs_comp_results()$g2_name) })
  output$ngs_comp_bar_plot <- renderPlot({ req(ngs_comp_results()); create_comp_barplot(ngs_comp_results()$summary) })
  
  observeEvent(input$ngs_comp_diff_click, {
    showModal(modalDialog(
      title = "NGS differential methylation (full size)",
      modal_fit_css,
      div(class = "panda-modal-plot",
          plotOutput("ngs_comp_diff_modal", width = "100%", height = "76vh")),
      size = "l", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$ngs_comp_diff_modal <- renderPlot({
    req(ngs_comp_results())
    create_diff_plot(ngs_comp_results()$site_table,
                     ngs_comp_results()$g1_name,
                     ngs_comp_results()$g2_name)
  }, res = 110)
  
  observeEvent(input$ngs_comp_bar_click, {
    showModal(modalDialog(
      title = "NGS group methylation comparison (full size)",
      modal_fit_css,
      div(class = "panda-modal-plot",
          plotOutput("ngs_comp_bar_modal", width = "100%", height = "76vh")),
      size = "l", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$ngs_comp_bar_modal <- renderPlot({
    req(ngs_comp_results())
    create_comp_barplot(ngs_comp_results()$summary)
  }, res = 110)
  output$ngs_comp_stat_txt <- renderPrint({ req(ngs_comp_results()); cat("P-value:", ngs_comp_results()$u_test_p) })
  
  output$ngs_comp_site_table <- renderDT({ 
    req(ngs_comp_results())
    comp <- ngs_comp_results()
    df <- comp$site_table
    colnames(df) <- c("Position", paste(comp$g1_name, "Meth"), paste(comp$g1_name, "Total"), paste(comp$g1_name, "Pct"), paste(comp$g2_name, "Meth"), paste(comp$g2_name, "Total"), paste(comp$g2_name, "Pct"), "P_Value", "FDR")
    datatable(df, rownames=FALSE) %>% formatRound(c(4, 7), 1) %>% formatSignif(c(8, 9), 3) 
  })
  
  # ==============================================================================
  # 4. DOWNLOAD HANDLERS & DYNAMIC UI
  # ==============================================================================
  
  output$ui_dl_sanger_single <- renderUI({
    if (is.null(sanger_single_res()) || is.null(sanger_single_het())) {
      actionButton("dummy_dl_ss", "Download Results (.zip)", icon = icon("download"), class = "w-100 btn-secondary disabled")
    } else {
      downloadButton("dl_sanger_single", "Download Results (.zip)", class = "w-100 btn-info")
    }
  })
  
  output$ui_dl_sanger_multi <- renderUI({
    if (is.null(sanger_multi_batch())) {
      actionButton("dummy_dl_sm", "Download Batch Data (.zip)", icon = icon("download"), class = "w-100 btn-secondary disabled")
    } else {
      downloadButton("dl_sanger_multi", "Download Batch Data (.zip)", class = "w-100 btn-info")
    }
  })
  
  output$ui_dl_sanger_comp <- renderUI({
    if (is.null(sanger_multi_comp())) {
      actionButton("dummy_dl_sc", "Download Comparison (.zip)", icon = icon("download"), class = "w-100 btn-secondary disabled")
    } else {
      downloadButton("dl_sanger_comp", "Download Comparison (.zip)", class = "w-100 btn-outline-success")
    }
  })
  
  output$ui_dl_ngs_batch <- renderUI({
    if (is.null(ngs_batch_results())) {
      actionButton("dummy_dl_nb", "Download Batch Data (.zip)", icon = icon("download"), class = "w-100 btn-secondary disabled")
    } else {
      downloadButton("dl_ngs_batch", "Download Batch Data (.zip)", class = "w-100 btn-info")
    }
  })
  
  output$ui_dl_ngs_comp <- renderUI({
    if (is.null(ngs_comp_results())) {
      actionButton("dummy_dl_nc", "Download Comparison (.zip)", icon = icon("download"), class = "w-100 btn-secondary disabled")
    } else {
      downloadButton("dl_ngs_comp", "Download Comparison (.zip)", class = "w-100 btn-outline-success")
    }
  })
  
  output$dl_sanger_single <- downloadHandler(
    filename = function() { paste0("PANDA_Single_", Sys.Date(), ".zip") },
    content = function(file) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      fs <- c("Summary.csv", "LongData.csv", "Heterogeneity.csv", "CpG_Stats.csv")
      write.csv(sanger_single_res()$read_summary, "Summary.csv", row.names=F)
      write.csv(sanger_single_res()$long_data, "LongData.csv", row.names=F)
      write.csv(sanger_single_het()$scores, "Heterogeneity.csv", row.names=F)
      write.csv(calculate_quma_stats(sanger_single_res(), "Sanger")$cpg_table, "CpG_Stats.csv", row.names=F)
      
      p_lol <- create_lollipop_plot(sanger_single_res()$long_data)
      if(!is.null(p_lol)) { ggsave("Plot_Lollipop.pdf", p_lol, width=13, height=7.15); fs <- c(fs, "Plot_Lollipop.pdf") }
      
      p_asm <- create_asm_heatmap(sanger_single_het())
      if(!is.null(p_asm)) { ggsave("Plot_ASM.pdf", p_asm, width=8, height=8); fs <- c(fs, "Plot_ASM.pdf") }
      
      zip(file, fs)
    }
  )
  
  output$dl_sanger_multi <- downloadHandler(
    filename = function() { paste0("PANDA_Sanger_Batch_", Sys.Date(), ".zip") },
    content = function(file) {
      batch <- sanger_multi_batch()
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      write.csv(batch$summary, "Batch_Summary.csv", row.names=F)
      
      try({
        d_sum <- batch$summary %>% filter(is.finite(Meth))
        p_sum <- ggplot(d_sum, aes(x=Sample, y=Meth, fill=Sample)) +
          geom_col(colour="black", show.legend=FALSE) +
          scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20)) + theme_minimal(base_size=14)
        ggsave("Batch_Summary_Plot.png", p_sum, width=10, height=6, dpi=600)
      })
      
      plot_dir <- "Individual_Results"
      if(dir.exists(plot_dir)) unlink(plot_dir, recursive=T)
      dir.create(plot_dir)
      
      withProgress(message = 'Saving All Details...', value = 0, {
        nms <- names(batch$details)
        n <- length(nms)
        for(i in seq_along(nms)) {
          nm <- nms[i]
          safe_nm <- gsub("[^a-zA-Z0-9_]", "_", nm)
          
          tryCatch({
            res <- batch$details[[nm]]
            het <- batch$hets[[nm]]
            stats <- calculate_quma_stats(res, "Sanger")
            
            write.csv(res$read_summary, file.path(plot_dir, paste0(safe_nm, "_Reads.csv")), row.names=F)
            write.csv(stats$cpg_table, file.path(plot_dir, paste0(safe_nm, "_CpG_Stats.csv")), row.names=F)
            write.csv(het$scores, file.path(plot_dir, paste0(safe_nm, "_Heterogeneity.csv")), row.names=F)
            
            if(!is.null(res$alignments$pairwise)) {
              writeLines(unlist(res$alignments$pairwise), file.path(plot_dir, paste0(safe_nm, "_Alignments.txt")))
            }
            
            try({
              p1 <- create_lollipop_plot(res$long_data)
              if(!is.null(p1)) ggsave(file.path(plot_dir, paste0(safe_nm, "_Lollipop.png")), p1, width=8, height=6, dpi=600)
            })
            
            try({
              p2 <- create_asm_heatmap(het)
              if(!is.null(p2)) ggsave(file.path(plot_dir, paste0(safe_nm, "_ASM.png")), p2, width=6, height=6, dpi=600)
            })
            
          }, error = function(e) { message(paste("Skipping sample:", nm, e$message)) })
          
          incProgress(1/n)
        }
      })
      
      fs <- c("Batch_Summary.csv")
      if(file.exists("Batch_Summary_Plot.png")) fs <- c(fs, "Batch_Summary_Plot.png")
      fs <- c(fs, list.files(plot_dir, full.names=T, recursive=T))
      
      if(length(fs) > 0) zip(file, fs) else stop("No files generated.")
    }
  )
  
  output$dl_sanger_comp <- downloadHandler(
    filename = function() { paste0("PANDA_Sanger_Compare_", Sys.Date(), ".zip") },
    content = function(file) {
      comp <- sanger_multi_comp()
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      fs <- c()
      
      df <- comp$site_table
      colnames(df) <- c("Position", paste(comp$g1_name, "Meth"), paste(comp$g1_name, "Total"), paste(comp$g1_name, "Pct"), paste(comp$g2_name, "Meth"), paste(comp$g2_name, "Total"), paste(comp$g2_name, "Pct"), "P_Value", "FDR")
      write.csv(df, "Comparison_Sites.csv", row.names=F)
      fs <- c(fs, "Comparison_Sites.csv")
      
      p1 <- create_diff_plot(comp$site_table, comp$g1_name, comp$g2_name)
      if(!is.null(p1)) { ggsave("Comp_Diff_Plot.pdf", p1, width=8, height=6); fs <- c(fs, "Comp_Diff_Plot.pdf") }
      
      p2 <- create_comp_barplot(comp$summary)
      if(!is.null(p2)) {
        ggsave("Comp_Mean_Bar.pdf", p2, width=6.5, height=5.5)
        fs <- c(fs, "Comp_Mean_Bar.pdf")
      }
      
      zip(file, fs)
    }
  )
  
  output$dl_ngs_batch <- downloadHandler(
    filename = function() { paste0("PANDA_NGS_Batch_", Sys.Date(), ".zip") },
    content = function(file) {
      batch <- ngs_batch_results()
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      write.csv(batch$summary, "NGS_Batch_Summary.csv", row.names=F)
      try({
        d_sum <- batch$summary %>% filter(is.finite(Meth))
        p_sum <- ggplot(d_sum, aes(x=Sample, y=Meth, fill=Sample)) +
          geom_col(colour="black", show.legend=FALSE) +
          scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20)) +
          theme_minimal(base_size=14) + theme(axis.text.x=element_text(angle=45, hjust=1))
        ggsave("NGS_Batch_Summary_Plot.png", p_sum, width=10, height=6, dpi=600)
      })
      
      plot_dir <- "NGS_Individual_Results"
      if(dir.exists(plot_dir)) unlink(plot_dir, recursive=T)
      dir.create(plot_dir)
      
      withProgress(message = 'Saving NGS Details...', value = 0, {
        nms <- names(batch$details)
        n <- length(nms)
        for(i in seq_along(nms)) {
          nm <- nms[i]
          safe_nm <- gsub("[^a-zA-Z0-9_]", "_", nm)
          
          tryCatch({
            res <- batch$details[[nm]]
            het <- batch$hets[[nm]]
            stats <- calculate_quma_stats(res, "NGS")
            
            write.csv(res$read_summary, file.path(plot_dir, paste0(safe_nm, "_Reads.csv")), row.names=F)
            write.csv(stats$cpg_table, file.path(plot_dir, paste0(safe_nm, "_CpG_Stats.csv")), row.names=F)
            write.csv(het$scores, file.path(plot_dir, paste0(safe_nm, "_Heterogeneity.csv")), row.names=F)
            
            if(!is.null(res$alignments$pairwise)) {
              writeLines(unlist(res$alignments$pairwise), file.path(plot_dir, paste0(safe_nm, "_Alignments.txt")))
            }
            
            n_reads <- length(unique(res$long_data$ReadID))
            calc_h <- 4 + (n_reads * 0.15)
            pt_size <- if(n_reads > 30) 3 else 4
            txt_size <- if(n_reads > 30) 10 else 14
            
            p1 <- create_meth_histogram(res$long_data)
            if(!is.null(p1)) ggsave(file.path(plot_dir, paste0(safe_nm, "_Histogram.png")), p1, width=8, height=6, dpi=600)
            
            p2 <- create_abundance_heatmap(res$long_data)
            if(!is.null(p2)) {
              hm_height <- max(8, min(10, 6 + 0.0015 * dplyr::n_distinct(res$long_data$ReadID)))
              ggsave(file.path(plot_dir, paste0(safe_nm, "_Abundance_Heatmap.png")),
                     p2, width=13, height=hm_height, dpi=600)
            }
            
            p3 <- create_lollipop_plot(res$long_data, point_size = pt_size, text_size = txt_size)
            if(!is.null(p3)) ggsave(file.path(plot_dir, paste0(safe_nm, "_Lollipop.png")), p3, width=8, height=calc_h, limitsize = FALSE, dpi=600)
            
          }, error = function(e) { message(paste("Skipping NGS sample:", nm, e$message)) })
          
          incProgress(1/n)
        }
      })
      
      fs <- c("NGS_Batch_Summary.csv")
      if(file.exists("NGS_Batch_Summary_Plot.png")) fs <- c(fs, "NGS_Batch_Summary_Plot.png")
      fs <- c(fs, list.files(plot_dir, full.names=T, recursive=T))
      
      if(length(fs) > 0) zip(file, fs) else stop("No files generated.")
    }
  )
  
  output$dl_ngs_comp <- downloadHandler(
    filename = function() { paste0("PANDA_NGS_Compare_", Sys.Date(), ".zip") },
    content = function(file) {
      comp <- ngs_comp_results()
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      fs <- c()
      
      df <- comp$site_table
      colnames(df) <- c("Position", paste(comp$g1_name, "Meth"), paste(comp$g1_name, "Total"), paste(comp$g1_name, "Pct"), paste(comp$g2_name, "Meth"), paste(comp$g2_name, "Total"), paste(comp$g2_name, "Pct"), "P_Value", "FDR")
      write.csv(df, "NGS_Comparison_Sites.csv", row.names=F)
      fs <- c(fs, "NGS_Comparison_Sites.csv")
      
      p1 <- create_diff_plot(comp$site_table, comp$g1_name, comp$g2_name)
      if(!is.null(p1)) { ggsave("NGS_Comp_Diff.pdf", p1, width=8, height=6); fs <- c(fs, "NGS_Comp_Diff.pdf") }
      
      p2 <- create_comp_barplot(comp$summary)
      if(!is.null(p2)) {
        ggsave("NGS_Comp_Bar.pdf", p2, width=6.5, height=5.5)
        fs <- c(fs, "NGS_Comp_Bar.pdf")
      }
      
      zip(file, fs)
    }
  )
}

shinyApp(ui, server)
