library(shiny)
library(shinyjs)
library(Biostrings)
library(tidyverse)
library(DT)
library(bslib)
library(patchwork)
library(sangerseqR)

# ==============================================================================
# 0. GLOBAL SETTINGS
# ==============================================================================
# Maximum 300MB for large NGS FASTQ files
options(shiny.maxRequestSize = 300 * 1024^2)

# ==============================================================================
# 1. CORE LOGIC (Backend Functions)
# ==============================================================================

run_bisulfite_alignment <- function(genome_seq, reads_set, 
                                    min_identity = 90, 
                                    min_conversion = 95,
                                    return_alignments = TRUE) {
  
  if (is(genome_seq, "DNAStringSet")) genome_seq <- genome_seq[[1]]
  genome_seq <- DNAString(as.character(genome_seq))
  
  # ----------------------------------------------------------------------------
  # Alignment Scoring Parameters
  # ----------------------------------------------------------------------------
  sub_mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = FALSE)
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
        aln_F <- pairwiseAlignment(pattern = r_F_T, subject = gen_T, type = "local", 
                                   substitutionMatrix = sub_mat, gapOpening = gap_op, gapExtension = gap_ext)
        
        r_R <- as.character(reverseComplement(DNAString(r)))
        r_R_T <- chartr("C", "T", r_R)
        aln_R <- pairwiseAlignment(pattern = r_R_T, subject = gen_T, type = "local", 
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
        raw_idx <- 1
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
  
  for (i in seq_len(n_total_reads)) {
    if (i %% 500 == 0) gc()
    
    raw_read_seq <- reads_set[[i]]
    read_name <- names(reads_set)[i]
    
    read_F <- raw_read_seq; read_F_T <- chartr("C", "T", read_F)
    aln_F <- pairwiseAlignment(pattern = read_F_T, subject = genome_T, type = "local", 
                               substitutionMatrix = sub_mat, gapOpening = gap_op, gapExtension = gap_ext)
    
    read_R <- reverseComplement(raw_read_seq); read_R_T <- chartr("C", "T", read_R)
    aln_R <- pairwiseAlignment(pattern = read_R_T, subject = genome_T, type = "local", 
                               substitutionMatrix = sub_mat, gapOpening = gap_op, gapExtension = gap_ext)
    
    if (score(aln_F) >= score(aln_R)) {
      best_aln <- aln_F; final_read_seq <- read_F; strand <- "Forward"
    } else {
      best_aln <- aln_R; final_read_seq <- read_R; strand <- "Reverse"
    }
    
    aln_pat_converted <- as.character(pattern(best_aln))
    raw_bases <- strsplit(as.character(final_read_seq), "")[[1]]
    aln_template <- strsplit(aln_pat_converted, "")[[1]]
    
    reconstructed_pat_chars <- character(length(aln_template))
    raw_idx <- 1
    
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
      
      orig_g_base <- as.character(subseq(genome_seq_char, curr_g_pos, curr_g_pos))
      
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

calculate_heterogeneity <- function(res_obj) {
  df_long <- res_obj$long_data
  cpg_sites <- res_obj$genome_info$cpg_pos
  
  if(is.null(df_long) || nrow(df_long) == 0) {
    return(list(scores = data.frame(Metric = c("PDR (Discordance)", "Epipolymorphism", "qFDRP"), Value = c(0, 0, 0)), meth_mat = matrix(NA), clusters = NULL))
  }
  
  meth_mat <- df_long %>% 
    select(ReadID, Position, Methylation) %>% 
    distinct(ReadID, Position, .keep_all = TRUE) %>% 
    pivot_wider(names_from = Position, values_from = Methylation) %>% 
    column_to_rownames("ReadID")
  
  miss <- setdiff(as.character(cpg_sites), colnames(meth_mat))
  if(length(miss) > 0) for(c in miss) meth_mat[[c]] <- NA
  meth_mat <- meth_mat[, as.character(cpg_sites), drop=FALSE] 
  
  calc_pdr <- function(row) { vals <- na.omit(row); if(length(vals) < 4) return(NA); return(as.integer(length(unique(vals)) > 1)) }
  pdr_score <- mean(apply(meth_mat, 1, calc_pdr), na.rm=TRUE) * 100
  if(is.nan(pdr_score)) pdr_score <- 0
  
  epialleles <- character()
  if(length(cpg_sites) >= 4) {
    for(i in 1:(length(cpg_sites)-3)) {
      window_mat <- meth_mat[, i:(i+3), drop=FALSE]
      patterns <- apply(window_mat, 1, function(x) { if(any(is.na(x))) return(NA); paste(x, collapse="") })
      epialleles <- c(epialleles, na.omit(patterns))
    }
  }
  if(length(epialleles) > 0) { probs <- table(epialleles)/length(epialleles); epipoly <- 1 - sum(probs^2)
  } else { epipoly <- 0 }
  
  cluster_mat <- meth_mat; cluster_mat[is.na(cluster_mat)] <- 0.5
  asm_cluster <- rep(1, nrow(cluster_mat))
  if(nrow(cluster_mat) >= 2) {
    try({ 
      set.seed(11) # Fixed seed for deterministic clustering
      asm_cluster <- kmeans(cluster_mat, centers = 2)$cluster 
    }, silent=TRUE)
  }
  
  used_mat <- meth_mat
  if(nrow(used_mat) > 100) {
    set.seed(11) # Fixed seed for deterministic subsampling
    used_mat <- used_mat[sample(nrow(used_mat), 100), ] 
  }
  qfdrp_sum <- 0; n_pairs <- 0; reads <- rownames(used_mat)
  
  if(length(reads) > 1) {
    for(i in 1:(length(reads)-1)) {
      for(j in (i+1):length(reads)) {
        r1 <- used_mat[i, ]; r2 <- used_mat[j, ]
        valid <- !is.na(r1) & !is.na(r2); n_ov <- sum(valid)
        if(n_ov > 0) { qfdrp_sum <- qfdrp_sum + (sum(r1[valid]!=r2[valid])/n_ov); n_pairs <- n_pairs + 1 }
      }
    }
  }
  qfdrp_score <- if(n_pairs > 0) (qfdrp_sum / n_pairs) else 0
  
  list(scores = data.frame(Metric = c("PDR (Discordance)", "Epipolymorphism", "qFDRP"), Value = c(pdr_score, epipoly, qfdrp_score)), meth_mat = meth_mat, clusters = asm_cluster)
}

calculate_quma_stats <- function(res_obj, mode = "Sanger") {
  df_long <- res_obj$long_data; df_summary <- res_obj$read_summary %>% filter(!str_detect(Pattern, "^excluded"))
  cpg_sites <- res_obj$genome_info$cpg_pos
  if(nrow(df_summary) == 0) return(list(overall=0, sd_cpg=0, se_cpg=0, sd_seq=0, se_seq=0, cpg_table=data.frame()))
  
  if (mode == "NGS") {
    df_long <- df_long %>% mutate(Count = as.integer(str_extract(ReadID, "(?<=Count)\\d+")))
    overall_meth <- sum(df_long$Methylation * df_long$Count) / sum(df_long$Count) * 100
    cpg_stats <- df_long %>% group_by(Position) %>% summarise(Num_Meth=sum(Methylation*Count), Num_Total=sum(Count), Meth_Pct=sum(Methylation*Count)/sum(Count)*100)
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
  
  p_vals <- numeric(nrow(full_site_table))
  for(i in 1:nrow(full_site_table)) {
    m1 <- full_site_table$Meth_1[i]; u1 <- full_site_table$Total_1[i] - m1
    m2 <- full_site_table$Meth_2[i]; u2 <- full_site_table$Total_2[i] - m2
    if(full_site_table$Total_1[i] == 0 && full_site_table$Total_2[i] == 0) { p_vals[i] <- NA } 
    else { mat <- matrix(c(m1, u1, m2, u2), nrow = 2, byrow = TRUE); p_vals[i] <- fisher.test(mat)$p.value }
  }
  full_site_table$P_Value <- p_vals
  
  # Benjamini-Hochberg FDR correction
  full_site_table$FDR <- p.adjust(p_vals, method = "BH")
  
  get_read_pcts <- function(res_list) {
    unlist(lapply(res_list, function(r) {
      d <- r$read_summary %>% filter(!str_detect(Pattern, "^excluded"))
      if(nrow(d)==0) return(numeric(0)); return(d$Meth_Pct) 
    }))
  }
  reads_A <- get_read_pcts(res_list_A); reads_B <- get_read_pcts(res_list_B)
  u_test_p <- NA
  if(length(reads_A) > 0 && length(reads_B) > 0) {
    reads_A <- na.omit(reads_A); reads_B <- na.omit(reads_B)
    if(length(reads_A) > 0 && length(reads_B) > 0) { 
      u_test <- wilcox.test(reads_A, reads_B, exact=FALSE)
      u_test_p <- u_test$p.value 
    }
  }
  
  sum_df <- data.frame(Group = c(g1_name, g2_name), Mean = c(mean(site_A$Pct, na.rm=T), mean(site_B$Pct, na.rm=T)))
  
  return(list(site_table = full_site_table, combined_long = combined_long, u_test_p = u_test_p, summary = sum_df, g1_name=g1_name, g2_name=g2_name))
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
  
  df_summary <- long_data %>%
    group_by(ReadID) %>%
    summarise(
      Mean_Meth = mean(Methylation, na.rm = TRUE),
      Count = if(any(grepl("Count", ReadID))) as.integer(str_extract(first(ReadID), "(?<=Count)\\d+")) else 1
    ) %>%
    arrange(desc(Mean_Meth)) %>%
    mutate(
      ymax = cumsum(Count),
      ymin = lag(ymax, default = 0)
    )
  
  plot_data <- long_data %>%
    inner_join(df_summary, by = "ReadID")
  
  ggplot(plot_data) +
    geom_rect(aes(xmin = Position - 2.5, xmax = Position + 2.5, 
                  ymin = ymin, ymax = ymax, 
                  fill = factor(Methylation)), color = NA) +
    scale_fill_manual(values = c("0" = "lightblue", "1" = "firebrick"), 
                      labels = c("Unmethylated", "Methylated"), 
                      name = "Status") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    labs(
      title = "Abundance Heatmap (Weighted by Read Count)",
      subtitle = "Vertical thickness of each segment represents allele read count. X-axis shows physical distance.",
      x = "CpG Position (bp)",
      y = "Cumulative Read Count"
    ) +
    theme(
      text = element_text(size=14),
      panel.grid = element_blank()
    )
}

create_lollipop_plot <- function(long_data, point_size = 4, text_size = 14) {
  if(is.null(long_data) || nrow(long_data) == 0) return(NULL)
  
  df <- long_data %>%
    mutate(CpG_Index = as.numeric(factor(Position)))
  
  read_stats <- df %>%
    group_by(ReadID) %>%
    summarise(
      Mean_Meth = mean(Methylation, na.rm = TRUE),
      Count = if(any(grepl("Count", ReadID))) as.integer(str_extract(first(ReadID), "(?<=Count)\\d+")) else 1
    ) %>%
    arrange(desc(Mean_Meth), desc(Count)) 
  
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
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
    ) +
    labs(x = "CpG Position (Sequential)", y = if(is_sanger) "Clones / Reads" else "Alleles (Sorted by Meth%)")
  
  if (!is_sanger) {
    p1 <- p1 + theme(axis.text.y = element_blank())
    p2 <- ggplot(read_stats, aes(x = Count, y = ReadID)) +
      geom_col(fill = "steelblue") +
      geom_text(aes(label = Count), hjust = -0.1, size = 3) +
      theme_void() + 
      theme(plot.margin = margin(l = 10, r = 20)) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
      labs(title = "Read Count")
    return(p1 + p2 + plot_layout(widths = c(3, 1)))
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
  ggplot(summary_df, aes(x=Group, y=Mean, fill=Group)) + 
    geom_bar(stat="identity", width=0.6, show.legend = FALSE) + 
    ylim(0,100) + 
    theme_minimal() +
    labs(y = "Mean Methylation (%)", title = "Global Methylation Difference")
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
  title = "PANDA", 
  theme = bs_theme(bootswatch = "minty"),
  header = shinyjs::useShinyjs(),
  
  nav_panel(title = "Introduction",
            card(
              h3("PANDA: Phased ANalysis of DNA Amplicons"),
              p(
                em("A unified platform for Phased Bisulfite Analysis and Heterogeneity Quantification (Sanger & NGS)"),
                br(), 
                "The source code, demo data, and detailed documentation are publicly available on ",
                a("https://github.com/kubo-azu/PANDA", 
                  href = "https://github.com/kubo-azu/PANDA", 
                  target = "_blank"),
                ".",
                br(),
                "Please cite this paper if you use this app in your research:",
                br(),
                "Kubota, A., Kobayashi, H., & Tajima, A. (2026). PANDA: Read-Level Phased Analysis of DNA Amplicons for Methylation Studies. bioRxiv. ",
                a("https://doi.org/10.64898/2026.04.01.715790", 
                  href = "https://doi.org/10.64898/2026.04.01.715790", 
                  target = "_blank"),
                "."
              ),
              
              div(class="alert alert-info",
                  h4(icon("flask"), " Quick Start with Demo Data"),
                  p("Download our comprehensive demo dataset (.zip) containing Sanger and NGS samples."),
                  downloadButton("dl_demo_data", "Download Demo Data (PANDA_Demo.zip)", class="btn-primary w-100")
              ),
              
              h4(icon("layer-group"), " Overview"),
              tags$ul(
                tags$li(strong("Sanger Mode:"), " Analysis for clonal sequencing (e.g. TA cloning). Supports single-sample detail view and multi-sample group comparison."),
                tags$li(strong("NGS Mode:"), " Optimized for high-depth amplicon sequencing. Includes automatic dereplication and statistical comparison."),
                tags$li(strong("Multi-Target & Incremental Upload:"), " Reference supports Multi-FASTA. Files can be accumulated in the queue without resetting.")
              ),
              
              h4(icon("list-check"), " How to Use (Step-by-Step)"),
              accordion(
                open = "Step 1: Data Preparation",
                
                accordion_panel("Step 1: Data Preparation",
                                p("The quality of analysis depends on a correct Reference Genome."),
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
                                  tags$li("Save as .fasta or .txt. (Multi-FASTA supported).")
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
                                              tags$li(strong("Heterogeneity & ASM:"), " Contains the 'ASM Heatmap' (clustering) and heterogeneity scores (PDR, qFDRP)."),
                                              tags$li(strong("Sequence Info & Alignments:"), " Quality control metrics and pairwise text alignments to check mismatches directly on the sequence.")
                                            ),
                                            hr(),
                                            
                                            h5(icon("server"), "2. NGS Mode"),
                                            h6("High-depth amplicon analysis optimized for ASM/Imprinting detection."),
                                            tags$ul(
                                              tags$li(strong("Top Sequences (Lollipop):"), " Shows the most abundant unique alleles. CpG sites are ", strong("equally spaced (QUMA-style)"), " for high visibility, combined with a bar chart of read counts. Ideal for presentations."),
                                              tags$li(strong("ASM Profile (Main):"),
                                                      tags$ul(
                                                        tags$li(strong("Heterogeneity & Distribution:"), " Advanced metrics (PDR, Epipoly) and a histogram to check bimodality."),
                                                        tags$li(strong("Abundance Heatmap (Weighted):"), " Visualizes allele proportions. Unlike the Lollipop plot, the X-axis reflects the ", strong("actual genomic distance (bp)"), ", allowing you to assess physical read coverage and drop-outs.")
                                                      )
                                              ),
                                              tags$li(strong("Statistics & Alignments:"), " Detailed tables for CpG rates. The Alignments tab shows raw read sequences, perfect for confirming Motif Filter targets.")
                                            ),
                                            hr(),
                                            
                                            h5(icon("scale-balanced"), "3. Group Comparison"),
                                            p("Statistical comparison (Difference Plot, Mann-Whitney U, Fisher's Exact) between two biological groups across Sanger or NGS batches. ", 
                                              strong("P-values for single-CpG comparisons are adjusted for multiple testing using the Benjamini-Hochberg (FDR) method."))
                                  ),
                                  
                                  nav_panel("Heterogeneity Metrics",
                                            p("PANDA calculates advanced metrics to quantify disorder in methylation patterns (based on Scherer et al., Nucleic Acids Res, 2020)."),
                                            hr(),
                                            
                                            h5(strong("PDR (Proportion of Discordant Reads):")),
                                            p("The percentage of reads that are locally discordant (i.e., contain both methylated and unmethylated CpGs). High PDR indicates stochastic or disordered methylation turnover."),
                                            hr(),
                                            
                                            h5(strong("Epipolymorphism:")),
                                            p("Calculated based on the probability of observing specific epialleles (patterns of 4 consecutive CpGs). High values indicate a diverse population of methylation patterns."),
                                            hr(),
                                            
                                            h5(strong("qFDRP (quantitative Fraction of Discordant Read Pairs):")),
                                            p("Measures the dissimilarity between randomly chosen pairs of reads. A high qFDRP indicates that the population contains diverse methylation states (high heterogeneity) rather than a uniform pattern.")
                                  )
                                )
                ),
                
                accordion_panel("Technical Note: Alignment & Processing",
                                
                                h6("1. Alignment & Quality Control"),
                                p("PANDA processes data completely within R, eliminating the need for external aligners."),
                                tags$ul(
                                  tags$li(strong("Alignment Strategy:"), " Employs ", strong("Local Pairwise Alignment"), " (Smith-Waterman algorithm via the ", code("Biostrings"), " package) on ", em("in silico"), " C-to-T converted sequences. ",
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
                    sliderInput("ngs_top_n", "Top N Alleles", 10, 100, 30),
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
                nav_panel("Batch Summary", plotOutput("ngs_batch_qc_plot", height="400px"), DTOutput("ngs_batch_table")),
                
                nav_panel("Single Sample Analysis", 
                          uiOutput("ngs_detail_selector"), 
                          navset_card_tab(
                            nav_panel("ASM Profile", 
                                      layout_columns(col_widths=c(12, 12),
                                                     card(card_header("1. Heterogeneity Metrics (PDR, Epipoly, qFDRP)"), DTOutput("ngs_detail_het")),
                                                     card(card_header("2. Methylation Distribution (Histogram)"), plotOutput("ngs_dist_plot", height="350px")),
                                                     card(card_header("3. Abundance Heatmap (Weighted)"), plotOutput("ngs_abund_heatmap", height="600px"))
                                      )
                            ),
                            nav_panel("Top Sequences", 
                                      h5("Top Abundant Alleles (Ranked)"),
                                      plotOutput("ngs_lollipop_plot", height="550px")
                            ),
                            nav_panel("Statistics", 
                                      div(style = "overflow-y: auto; max-height: 800px;",
                                          h5("CpG Methylation Rates"), DTOutput("ngs_stats_table"),
                                          br(), hr(), br(),
                                          h5("Read Quality Control"), DTOutput("ngs_read_qc"),
                                          br(), br()
                                      )
                            ),
                            nav_panel("Alignments", 
                                      verbatimTextOutput("ngs_alignments_text")
                            )
                          )
                ),
                
                nav_panel("Group Comparison", h5("Differential Methylation"), layout_columns(col_widths=c(6,6), plotOutput("ngs_comp_diff_plot"), plotOutput("ngs_comp_bar_plot")), hr(), h5("Statistics"), verbatimTextOutput("ngs_comp_stat_txt"), DTOutput("ngs_comp_site_table"))
              )
            )
  )
)

# ==============================================================================
# 3. SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  output$dl_demo_data <- downloadHandler(
    filename = function() { "PANDA_Demo.zip" },
    content = function(file) {
      if(file.exists("PANDA_Demo.zip")) {
        file.copy("PANDA_Demo.zip", file)
      } else {
        showNotification("Demo data file not found in directory.", type="error")
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
    new_df <- input$sanger_multi_files; old_df <- stored_sanger(); stored_sanger(bind_rows(old_df, new_df))
  })
  
  # Clear Sanger Files
  observeEvent(input$clear_sanger_files, { 
    shinyjs::reset("sanger_multi_files")
    stored_sanger(data.frame()) 
    sanger_multi_batch(NULL)
    sanger_multi_comp(NULL)
  })
  output$sanger_file_count_txt <- renderText({ n <- nrow(stored_sanger()); if(n==0) "No files selected." else paste(n, "files stored.") })
  
  observeEvent(input$ngs_files_merged, {
    new_df <- input$ngs_files_merged
    if(!is.null(new_df)) {
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
      sanger_genomes_list(readDNAStringSet(input$sanger_genome$datapath))
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
          reads <- readDNAStringSet(input$sanger_single_fasta$datapath) 
        } else { 
          req(input$sanger_single_ab1); 
          reads <- process_ab1_files(input$sanger_single_ab1$datapath, input$sanger_single_ab1$name, input$ab1_trim_start, input$ab1_trim_end) 
        }
        incProgress(0.6, detail = "Aligning..."); res <- run_bisulfite_alignment(genome, reads, input$sanger_ident, input$sanger_conv, return_alignments = TRUE) 
        if(isTRUE(res$was_flipped)) showNotification("Reference automatically flipped to Reverse Complement.", type="warning", duration=10)
        
        sanger_single_res(res); sanger_single_het(calculate_heterogeneity(res)); incProgress(1, detail = "Done!")
        
      }, error = function(e) showNotification(e$message, type="error"))
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
    req(get_sanger_target_seq()); files <- stored_sanger(); if(nrow(files)==0) { showNotification("No files selected.", type="error"); return() }
    genome <- get_sanger_target_seq(); res_list <- list(); het_list <- list(); sum_list <- list()
    
    target_motifs <- character()
    if (!is.null(input$sanger_motif_file)) {
      tryCatch({
        target_motifs <- c(target_motifs, readLines(input$sanger_motif_file$datapath))
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
      showNotification("No valid results. Check file format or filter criteria.", type="error") 
    }
  })
  
  output$sanger_batch_qc <- renderPlot({ req(sanger_multi_batch()); ggplot(sanger_multi_batch()$summary, aes(x=Sample, y=Meth, fill=Sample)) + geom_bar(stat="identity", col="black", show.legend=FALSE) + ylim(0,100) + theme_minimal() })
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
      showNotification("Error: Please run Batch Analysis first.", type = "error")
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
      ngs_genomes_list(readDNAStringSet(input$ngs_genome$datapath))
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
    req(ngs_genomes_list(), input$ngs_target_select); files <- stored_ngs(); if(nrow(files)==0) { showNotification("No files.", type="error"); return() }
    
    all_genomes <- ngs_genomes_list()
    if(is.null(names(all_genomes))) names(all_genomes) <- paste0("Target_", seq_along(all_genomes))
    target_genome <- all_genomes[[input$ngs_target_select]]
    
    target_motifs <- character()
    if (!is.null(input$ngs_motif_file)) {
      tryCatch({
        target_motifs <- c(target_motifs, readLines(input$ngs_motif_file$datapath))
      }, error = function(e) showNotification("Error reading motif file", type="warning"))
    }
    if (input$ngs_motif_text != "") {
      target_motifs <- c(target_motifs, strsplit(input$ngs_motif_text, ",")[[1]])
    }
    target_motifs <- trimws(target_motifs)
    target_motifs <- target_motifs[target_motifs != ""]
    
    res_list <- list(); het_list <- list(); sum_list <- list()
    
    withProgress(message = 'NGS Batch Analysis...', value = 0, {
      n <- nrow(files)
      for(i in 1:n) {
        incProgress(1/n, detail = files$name[i])
        
        tryCatch({
          read_seqs <- function(fp, orig_name) {
            if (is.na(fp) || is.na(orig_name)) return(character())
            con <- if (grepl("\\.gz$", orig_name, ignore.case = TRUE)) gzfile(fp, "rt") else file(fp, "rt")
            lines <- readLines(con)
            close(con)
            
            if (input$ngs_type == "fastq") {
              if(length(lines) >= 2) return(lines[seq(2, length(lines), by = 4)]) else return(character())
            } else {
              return(lines[!grepl("^>", lines)])
            }
          }
          
          reads_char <- read_seqs(files$datapath[i], files$orig_name_R1[i])
          is_unmerged_paired <- (input$ngs_read_mode == "unmerged" && !is.na(files$datapath_R2[i]))
          
          if (is_unmerged_paired) {
            reads_char_r2 <- read_seqs(files$datapath_R2[i], files$orig_name_R2[i])
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
                reads_set <- DNAStringSet()
              }
            } else {
              reads_set <- DNAStringSet()
            }
          } else {
            if (!is.na(files$datapath_R2[i])) {
              reads_char_r2 <- read_seqs(files$datapath_R2[i], files$orig_name_R2[i])
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
              top_seqs <- head(counts, input$ngs_top_n)
              
              reads_set <- DNAStringSet(names(top_seqs))
              names(reads_set) <- paste0("Rank", 1:length(top_seqs), "_Count", as.integer(top_seqs))
            } else {
              reads_set <- DNAStringSet()
            }
          }
          
          if (length(reads_set) > 0) {
            res <- run_bisulfite_alignment(target_genome, reads_set, input$ngs_ident, input$ngs_conv, return_alignments = TRUE) 
            
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
              
              nm <- files$name[i]
              res_list[[nm]] <- res
              het_list[[nm]] <- calculate_heterogeneity(res)
              stats <- calculate_quma_stats(res, "NGS")
              
              real_total_reads <- sum(res$read_summary$Count)
              real_used_reads <- sum(res$read_summary$Count[res$read_summary$Pattern == "Passed" | !grepl("^excluded", res$read_summary$Pattern)])
              
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
          showNotification(paste("Error in", files$name[i], ":", e$message), type="error") 
        })
      }
    })
    
    if(length(res_list)>0) { 
      ngs_batch_results(list(details=res_list, hets=het_list, summary=bind_rows(sum_list), genome=target_genome)) 
      if(length(target_motifs) > 0) showNotification(paste("Filtered by", length(target_motifs), "motifs."), type="message")
    } else { 
      showNotification("No valid results. Check file format or filter criteria.", type="error") 
    }
  })
  
  output$ngs_batch_qc_plot <- renderPlot({ req(ngs_batch_results()); if(nrow(ngs_batch_results()$summary)>0) ggplot(ngs_batch_results()$summary, aes(x=Sample, y=Meth, fill=Sample)) + geom_bar(stat="identity", col="black", show.legend=FALSE) + ylim(0,100) + theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1)) })
  output$ngs_batch_table <- renderDT({ req(ngs_batch_results()); datatable(ngs_batch_results()$summary, rownames=FALSE) })
  output$ngs_detail_selector <- renderUI({ req(ngs_batch_results()); selectInput("ngs_view_sample", "Select Sample:", choices=names(ngs_batch_results()$details)) })
  get_ngs_sel <- reactive({ req(ngs_batch_results(), input$ngs_view_sample); list(res = ngs_batch_results()$details[[input$ngs_view_sample]], het = ngs_batch_results()$hets[[input$ngs_view_sample]]) })
  
  output$ngs_dist_plot <- renderPlot({ req(get_ngs_sel()); create_meth_histogram(get_ngs_sel()$res$long_data) })
  output$ngs_abund_heatmap <- renderPlot({ req(get_ngs_sel()); create_abundance_heatmap(get_ngs_sel()$res$long_data) })
  output$ngs_lollipop_plot <- renderPlot({ req(get_ngs_sel()); create_lollipop_plot(get_ngs_sel()$res$long_data) })
  
  output$ngs_detail_het <- renderDT({ datatable(get_ngs_sel()$het$scores, options=list(dom='t'), rownames=FALSE) %>% formatRound(2, 3) })
  
  output$ngs_stats_table <- renderDT({ datatable(calculate_quma_stats(get_ngs_sel()$res, "NGS")$cpg_table, rownames=FALSE) %>% formatRound("Meth_Pct", 1) })
  output$ngs_read_qc <- renderDT({ datatable(get_ngs_sel()$res$read_summary, options=list(scrollX=T), rownames=FALSE) })
  
  output$ngs_alignments_text <- renderPrint({ req(get_ngs_sel()); cat(paste(get_ngs_sel()$res$alignments$pairwise, collapse="\n----------------\n")) })
  
  output$ngs_group_ui <- renderUI({ req(ngs_batch_results()); tagList(selectInput("ngs_grp1", "G1:", names(ngs_batch_results()$details), multiple=T), selectInput("ngs_grp2", "G2:", names(ngs_batch_results()$details), multiple=T)) })
  
  # Run NGS Comp
  observeEvent(input$run_ngs_comp, { 
    if (is.null(ngs_batch_results())) {
      showNotification("Error: Please run Batch Analysis first.", type = "error")
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
      if(!is.null(p_lol)) { ggsave("Plot_Lollipop.pdf", p_lol, width=10, height=5.5); fs <- c(fs, "Plot_Lollipop.pdf") }
      
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
        p_sum <- ggplot(batch$summary, aes(x=Sample, y=Meth, fill=Sample)) + 
          geom_bar(stat="identity", col="black", show.legend=FALSE) + ylim(0,100) + theme_minimal()
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
      if(!is.null(p2)) { ggsave("Comp_Mean_Bar.pdf", p2, width=6, height=6); fs <- c(fs, "Comp_Mean_Bar.pdf") }
      
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
        p_sum <- ggplot(batch$summary, aes(x=Sample, y=Meth, fill=Sample)) + 
          geom_bar(stat="identity", col="black", show.legend=FALSE) + ylim(0,100) + theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1))
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
            if(!is.null(p2)) ggsave(file.path(plot_dir, paste0(safe_nm, "_Abundance_Heatmap.png")), p2, width=8, height=6, dpi=600)
            
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
      if(!is.null(p2)) { ggsave("NGS_Comp_Bar.pdf", p2, width=6, height=6); fs <- c(fs, "NGS_Comp_Bar.pdf") }
      
      zip(file, fs)
    }
  )
}

shinyApp(ui, server)
