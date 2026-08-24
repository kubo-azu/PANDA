# PANDA alignment functions

# pwalign is required by the current implementation.  The deprecated
# Biostrings alignment API is intentionally not used as a fallback.
.panda_has_pwalign <- TRUE

# Alignment objects are cached within the current R process.  This is
# particularly important for GUI NGS analysis, where the Top-N display set
# overlaps with the full dereplicated set used for quantitative metrics.
# The cache changes no results; it only avoids repeating identical calls.
.panda_alignment_cache <- new.env(parent = emptyenv())
.panda_alignment_cache_limit <- 10000L

.panda_pairwiseAlignment <- function(...) {
  args <- list(...)
  pattern_key <- as.character(args$pattern)
  subject_key <- as.character(args$subject)
  type_key <- if (is.null(args$type)) "" else as.character(args$type)
  gap_open_key <- if (is.null(args$gapOpening)) "" else as.character(args$gapOpening)
  gap_ext_key <- if (is.null(args$gapExtension)) "" else as.character(args$gapExtension)
  key <- paste(
    type_key, gap_open_key, gap_ext_key,
    pattern_key, subject_key,
    sep = "\r"
  )
  
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

run_bisulfite_alignment <- function(genome_seq, reads_set, min_identity = 90,
                                    min_conversion = 95, return_alignments = TRUE, workers = 1L) {
  if (methods::is(genome_seq, "DNAStringSet")) 
    genome_seq <- genome_seq[[1]]
  genome_seq <- DNAString(as.character(genome_seq))
  sub_mat <- .panda_nucleotideSubstitutionMatrix(match = 1, 
                                                 mismatch = -3, baseOnly = FALSE)
  gap_op <- -10
  gap_ext <- -4
  was_flipped <- FALSE
  n_test <- min(10, length(reads_set))
  if (n_test > 0) {
    test_reads <- reads_set[1:n_test]
    test_conv <- function(test_gen) {
      test_gen_char <- as.character(test_gen)
      gen_T <- chartr("C", "T", test_gen_char)
      c_conv <- 0
      c_unconv <- 0
      for (i in seq_len(n_test)) {
        r <- as.character(test_reads[[i]])
        r_F <- r
        r_F_T <- chartr("C", "T", r_F)
        aln_F <- .panda_pairwiseAlignment(pattern = r_F_T, 
                                          subject = gen_T, type = "global", substitutionMatrix = sub_mat, 
                                          gapOpening = gap_op, gapExtension = gap_ext)
        r_R <- as.character(reverseComplement(DNAString(r)))
        r_R_T <- chartr("C", "T", r_R)
        aln_R <- .panda_pairwiseAlignment(pattern = r_R_T, 
                                          subject = gen_T, type = "global", substitutionMatrix = sub_mat, 
                                          gapOpening = gap_op, gapExtension = gap_ext)
        if (score(aln_F) >= score(aln_R)) {
          best_aln <- aln_F
          final_r <- r_F
        }
        else {
          best_aln <- aln_R
          final_r <- r_R
        }
        aln_pat_converted <- as.character(pattern(best_aln))
        raw_bases <- strsplit(final_r, "")[[1]]
        aln_template <- strsplit(aln_pat_converted, "")[[1]]
        reconstructed <- character(length(aln_template))
        raw_idx <- start(pattern(best_aln))
        for (j in seq_along(aln_template)) {
          if (aln_template[j] == "-") 
            reconstructed[j] <- "-"
          else {
            if (raw_idx <= length(raw_bases)) {
              reconstructed[j] <- raw_bases[raw_idx]
              raw_idx <- raw_idx + 1
            }
            else reconstructed[j] <- "N"
          }
        }
        aln_sub_str <- as.character(subject(best_aln))
        start_genome <- start(subject(best_aln))
        curr_g_pos <- start_genome - 1
        sub_chars <- strsplit(aln_sub_str, "")[[1]]
        for (k in seq_along(sub_chars)) {
          s_char <- sub_chars[k]
          p_char <- reconstructed[k]
          if (s_char != "-") 
            curr_g_pos <- curr_g_pos + 1
          if (s_char == "-" || p_char == "-") 
            next
          orig_g_base <- substring(test_gen_char, curr_g_pos, 
                                   curr_g_pos)
          if (orig_g_base == "C") {
            if (p_char == "T") 
              c_conv <- c_conv + 1
            else if (p_char == "C") 
              c_unconv <- c_unconv + 1
          }
        }
      }
      if ((c_conv + c_unconv) == 0) 
        return(0)
      else return((c_conv/(c_conv + c_unconv)) * 100)
    }
    conv_fwd <- test_conv(genome_seq)
    conv_rc <- test_conv(reverseComplement(genome_seq))
    if (conv_rc > 50 && conv_rc > (conv_fwd + 20)) {
      genome_seq <- reverseComplement(genome_seq)
      was_flipped <- TRUE
    }
  }
  genome_seq_char <- as.character(genome_seq)
  cpg_hits <- matchPattern("CG", genome_seq_char)
  cpg_sites <- start(cpg_hits)
  genome_T <- chartr("C", "T", genome_seq_char)
  res_ids <- character()
  res_strands <- character()
  res_pos <- integer()
  res_meth <- integer()
  sum_id <- character()
  sum_strand <- character()
  sum_mm <- integer()
  sum_gaps <- integer()
  sum_ident <- numeric()
  sum_meth_pct <- numeric()
  sum_conv_pct <- numeric()
  sum_pattern <- character()
  sum_meth_cpgs <- integer()
  pairwise_txt_list <- list()
  multi_align_seqs <- list()
  n_total_reads <- length(reads_set)
  n_excluded <- 0
  
  workers <- suppressWarnings(as.integer(workers))
  if (is.na(workers) || workers < 1L) workers <- 1L
  if (workers > 1L && .Platform$OS.type == "windows") {
    message("  Windows detected: using serial alignment for portability.")
    workers <- 1L
  }
  
  ## Precompute the two strand alignments in parallel.  The alignment
  ## algorithm and scoring are unchanged; only scheduling differs.
  align_one_read <- function(i) {
    raw_read_seq <- reads_set[[i]]
    read_F <- raw_read_seq
    read_F_T <- chartr("C", "T", read_F)
    aln_F <- .panda_pairwiseAlignment(
      pattern = read_F_T, subject = genome_T, type = "global",
      substitutionMatrix = sub_mat,
      gapOpening = gap_op, gapExtension = gap_ext
    )
    read_R <- reverseComplement(raw_read_seq)
    read_R_T <- chartr("C", "T", read_R)
    aln_R <- .panda_pairwiseAlignment(
      pattern = read_R_T, subject = genome_T, type = "global",
      substitutionMatrix = sub_mat,
      gapOpening = gap_op, gapExtension = gap_ext
    )
    list(read_F = read_F, read_R = read_R, aln_F = aln_F, aln_R = aln_R)
  }
  
  alignment_pairs <- if (workers > 1L && n_total_reads > 1L) {
    message(
      "  Parallel alignment: ", n_total_reads,
      " representatives on ", workers, " workers"
    )
    parallel::mclapply(
      seq_len(n_total_reads), align_one_read,
      mc.cores = workers, mc.preschedule = TRUE
    )
  } else {
    lapply(seq_len(n_total_reads), align_one_read)
  }
  
  for (i in seq_len(n_total_reads)) {
    if (i == 1L || i %% 100L == 0L || i == n_total_reads) {
      message(
        "  Aligning representative ", i, "/", n_total_reads
      )
    }
    if (i%%500 == 0)
      gc()
    alignment_pair <- alignment_pairs[[i]]
    raw_read_seq <- alignment_pair$read_F
    read_name <- names(reads_set)[i]
    read_F <- alignment_pair$read_F
    read_R <- alignment_pair$read_R
    aln_F <- alignment_pair$aln_F
    aln_R <- alignment_pair$aln_R
    if (score(aln_F) >= score(aln_R)) {
      best_aln <- aln_F
      final_read_seq <- read_F
      strand <- "Forward"
    }
    else {
      best_aln <- aln_R
      final_read_seq <- read_R
      strand <- "Reverse"
    }
    aln_pat_converted <- as.character(pattern(best_aln))
    raw_bases <- strsplit(as.character(final_read_seq), "")[[1]]
    aln_template <- strsplit(aln_pat_converted, "")[[1]]
    reconstructed_pat_chars <- character(length(aln_template))
    raw_idx <- start(pattern(best_aln))
    for (j in seq_along(aln_template)) {
      if (aln_template[j] == "-") {
        reconstructed_pat_chars[j] <- "-"
      }
      else {
        if (raw_idx <= length(raw_bases)) {
          reconstructed_pat_chars[j] <- raw_bases[raw_idx]
          raw_idx <- raw_idx + 1
        }
        else {
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
    identity_score <- (n_match/aln_len) * 100
    start_genome <- start(subject(best_aln))
    curr_g_pos <- start_genome - 1
    sub_chars <- strsplit(aln_sub_str, "")[[1]]
    pat_chars <- reconstructed_pat_chars
    tmp_meth <- integer()
    tmp_pos <- integer()
    conv_C_count <- 0
    unconv_C_count <- 0
    for (k in 1:aln_len) {
      s_char <- sub_chars[k]
      p_char <- pat_chars[k]
      if (s_char != "-") 
        curr_g_pos <- curr_g_pos + 1
      if (s_char == "-" || p_char == "-") 
        next
      ## genome_seq_char is a single character string; use substring()
      ## for one-based genomic coordinates.
      orig_g_base <- substring(genome_seq_char, curr_g_pos, curr_g_pos)
      if (orig_g_base == "C") {
        if (curr_g_pos %in% cpg_sites) {
          if (p_char == "C") {
            tmp_meth <- c(tmp_meth, 1L)
            tmp_pos <- c(tmp_pos, curr_g_pos)
          }
          else if (p_char == "T") {
            tmp_meth <- c(tmp_meth, 0L)
            tmp_pos <- c(tmp_pos, curr_g_pos)
          }
        }
        else {
          if (p_char == "C") 
            unconv_C_count <- unconv_C_count + 1
          else if (p_char == "T") 
            conv_C_count <- conv_C_count + 1
        }
      }
    }
    total_cph <- unconv_C_count + conv_C_count
    conv_rate <- if (total_cph > 0) 
      (conv_C_count/total_cph) * 100
    else 100
    exclusion_reason <- ""
    if (identity_score < min_identity) 
      exclusion_reason <- paste0("excluded (Id:", round(identity_score, 
                                                        1), "%)")
    else if (conv_rate < min_conversion) 
      exclusion_reason <- paste0("excluded (Conv:", round(conv_rate, 
                                                          1), "%)")
    if (exclusion_reason != "") 
      n_excluded <- n_excluded + 1
    sum_id <- c(sum_id, read_name)
    sum_strand <- c(sum_strand, strand)
    sum_mm <- c(sum_mm, aln_len - n_match)
    sum_gaps <- c(sum_gaps, total_gaps)
    sum_ident <- c(sum_ident, round(identity_score, 1))
    if (exclusion_reason == "") {
      m_pct <- if (length(tmp_meth) > 0) 
        round(mean(tmp_meth) * 100, 1)
      else NA
      sum_meth_pct <- c(sum_meth_pct, m_pct)
      sum_conv_pct <- c(sum_conv_pct, round(conv_rate, 
                                            1))
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
        pair_txt <- paste0("> ", read_name, "\nGen: ", 
                           aln_sub_str, "\nSeq: ", paste(reconstructed_pat_chars, 
                                                         collapse = ""), "\n")
        pairwise_txt_list[[length(pairwise_txt_list) + 
                             1]] <- pair_txt
      }
    }
    else {
      sum_meth_pct <- c(sum_meth_pct, NA)
      sum_conv_pct <- c(sum_conv_pct, round(conv_rate, 
                                            1))
      sum_pattern <- c(sum_pattern, exclusion_reason)
      sum_meth_cpgs <- c(sum_meth_cpgs, 0)
    }
  }
  if (length(sum_id) == 0) 
    return(NULL)
  read_summary_df <- data.frame(ReadID = sum_id, Strand = sum_strand, 
                                Mismatches = sum_mm, Gaps = sum_gaps, Identity_Pct = sum_ident, 
                                Meth_Pct = sum_meth_pct, Conv_Pct = sum_conv_pct, CpG_Count = sum_meth_cpgs, 
                                Pattern = sum_pattern, stringsAsFactors = FALSE)
  long_data_df <- data.frame(ReadID = res_ids, Strand = res_strands, 
                             Position = res_pos, Methylation = res_meth, stringsAsFactors = FALSE)
  return(list(long_data = long_data_df, read_summary = read_summary_df, 
              genome_info = list(len = nchar(genome_seq_char), n_cpg = length(cpg_sites), 
                                 cpg_pos = cpg_sites, seq = genome_seq_char), counts = list(total = n_total_reads, 
                                                                                            used = n_total_reads - n_excluded, excluded = n_excluded), 
              alignments = list(pairwise = pairwise_txt_list, multi = multi_align_seqs), 
              was_flipped = was_flipped))
}

process_ab1_files <- function(file_paths, file_names, trim_start = 20, 
                              trim_end = 20) {
  seq_list <- DNAStringSet()
  for (i in seq_along(file_paths)) {
    tryCatch({
      sanger <- readsangerseq(file_paths[i])
      seq <- primarySeq(sanger)
      if (length(seq) > (trim_start + trim_end)) 
        seq <- subseq(seq, start = trim_start + 1, end = length(seq) - 
                        trim_end)
      current_set <- DNAStringSet(seq)
      names(current_set) <- file_names[i]
      seq_list <- c(seq_list, current_set)
    }, error = function(e) warning(paste("Failed:", file_names[i])))
  }
  return(seq_list)
}

.panda_read_sequences <- function(filepath, filename, format = c("fastq", 
                                                                 "fasta")) {
  format <- match.arg(format)
  if (is.na(filepath) || is.na(filename) || !file.exists(filepath)) 
    return(character())
  con <- if (grepl("\\.gz$", filename, ignore.case = TRUE)) 
    gzfile(filepath, "rt")
  else file(filepath, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  if (!length(lines)) 
    return(character())
  if (format == "fastq") {
    if (length(lines)%%4L != 0L || any(!startsWith(lines[seq(1L, 
                                                             length(lines), by = 4L)], "@"))) 
      stop("Invalid FASTQ structure: expected 4 lines per record.")
    seqs <- lines[seq(2L, length(lines), by = 4L)]
    quals <- lines[seq(4L, length(lines), by = 4L)]
    if (any(nchar(seqs) != nchar(quals))) 
      stop("FASTQ sequence/quality length mismatch.")
  }
  else {
    hdr <- which(startsWith(lines, ">"))
    if (!length(hdr)) 
      stop("FASTA file contains no header.")
    ends <- c(hdr[-1L] - 1L, length(lines))
    seqs <- vapply(seq_along(hdr), function(ii) {
      x <- lines[(hdr[ii] + 1L):ends[ii]]
      paste0(x[!startsWith(x, ">")], collapse = "")
    }, character(1))
  }
  seqs <- toupper(gsub("\\s+", "", seqs))
  if (any(!grepl("^[ACGTN]+$", seqs))) 
    stop("Sequence contains unsupported bases.")
  seqs[nzchar(seqs)]
}

.panda_expand_dereplicated <- function(counts, top_n = Inf, paired = FALSE) {
  if (!length(counts)) 
    return(DNAStringSet())
  counts <- utils::head(sort(counts, decreasing = TRUE), top_n)
  if (!paired) {
    out <- DNAStringSet(names(counts))
    names(out) <- paste0("Rank", seq_along(counts), "_Count", 
                         as.integer(counts))
    return(out)
  }
  seqs <- character()
  ids <- character()
  for (ii in seq_along(counts)) {
    parts <- strsplit(names(counts)[ii], "---PAIR---", fixed = TRUE)[[1L]]
    base <- paste0("Rank", ii, "_Count", as.integer(counts[ii]))
    seqs <- c(seqs, parts[1L], parts[2L])
    ids <- c(ids, paste0(base, "_R1"), paste0(base, "_R2"))
  }
  out <- DNAStringSet(seqs)
  names(out) <- ids
  out
}

.panda_finalize_ngs_result <- function(res, is_unmerged = FALSE) {
  if (is.null(res)) 
    return(NULL)
  if (is_unmerged) {
    res$read_summary <- res$read_summary %>% mutate(BaseID = sub("_R[12]$", 
                                                                 "", ReadID)) %>% group_by(BaseID) %>% summarise(Strand = paste(unique(Strand), 
                                                                                                                                collapse = "/"), Mismatches = sum(Mismatches), Gaps = sum(Gaps), 
                                                                                                                 Identity_Pct = mean(Identity_Pct), Meth_Pct = mean(Meth_Pct, 
                                                                                                                                                                    na.rm = TRUE), Conv_Pct = mean(Conv_Pct), CpG_Count = sum(CpG_Count), 
                                                                                                                 Pattern = if (any(grepl("excluded", Pattern))) 
                                                                                                                   "excluded (Pair Failed QC)"
                                                                                                                 else "Passed", .groups = "drop") %>% rename(ReadID = BaseID)
    res$long_data <- res$long_data %>% mutate(ReadID = sub("_R[12]$", 
                                                           "", ReadID))
    excluded_ids <- res$read_summary$ReadID[grepl("excluded", 
                                                  res$read_summary$Pattern)]
    res$long_data <- res$long_data %>% filter(!ReadID %in% 
                                                excluded_ids)
    res$counts$total <- res$counts$total/2
  }
  res$long_data <- res$long_data %>% mutate(Count = as.integer(str_extract(ReadID, 
                                                                           "(?<=Count)\\d+")), Count = replace_na(Count, 1L))
  res$read_summary <- res$read_summary %>% mutate(Count = as.integer(str_extract(ReadID, 
                                                                                 "(?<=Count)\\d+")), Count = replace_na(Count, 1L))
  res
}
