# PANDA heterogeneity and statistical functions
# Extracted from app.R during CLI refactoring.

# Exact abundance-weighted qFDRP calculation using coverage-pattern
# aggregation. Count is an observed read-abundance weight, not an inferred
# number of original biological molecules. The calculation is exactly
# equivalent to expanding each dereplicated NGS variant Count times, including
# zero-distance pairs formed by two reads of the same variant.
.panda_weighted_median <- function(values, weights) {
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(keep)) return(NA_real_)
  values <- values[keep]
  weights <- weights[keep]
  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord]
  values[which(cumsum(weights) >= sum(weights) / 2)[[1L]]]
}

.panda_weighted_qfdrp <- function(meth_mat, weights, min_shared_cpg = 4L) {
  n <- nrow(meth_mat)
  min_shared_cpg <- suppressWarnings(as.integer(min_shared_cpg[[1L]]))
  if (is.na(min_shared_cpg) || min_shared_cpg < 1L) {
    stop("min_shared_cpg must be a positive integer.", call. = FALSE)
  }
  if (n < 1L || ncol(meth_mat) < 1L) {
    return(list(
      value = NA_real_, status = "no_eligible_pairs",
      n_eligible_variant_pairs = 0, weighted_eligible_pairs = 0,
      median_shared_cpg = NA_real_, shared_cpg_distribution = data.frame()
    ))
  }

  weights <- as.numeric(weights)
  weights[is.na(weights) | weights < 1] <- 1
  masks <- apply(!is.na(meth_mat), 1L, paste0, collapse = "")
  groups <- split(seq_len(n), masks)
  q_sum <- 0
  q_den <- 0
  n_variant_pairs <- 0
  shared_values <- numeric()
  shared_weights <- numeric()

  add_shared_diagnostic <- function(shared_cpg, pair_weight) {
    shared_values <<- c(shared_values, shared_cpg)
    shared_weights <<- c(shared_weights, pair_weight)
  }

  for (aa in seq_along(groups)) {
    ia <- groups[[aa]]
    xa <- meth_mat[ia, , drop = FALSE]
    wa <- weights[ia]
    total_a <- sum(wa)
    shared_self <- !is.na(meth_mat[ia[[1L]], ])
    k_self <- sum(shared_self)

    # All pairs represented by one coverage-mask group are eligible together.
    # choose(total_a, 2) includes pairs of reads belonging to the same exact
    # variant; their Hamming distance is zero and therefore they contribute to
    # the denominator but not to the numerator.
    if (k_self >= min_shared_cpg && total_a >= 2) {
      pair_weight <- total_a * (total_a - 1) / 2
      meth_a <- colSums(xa[, shared_self, drop = FALSE] * wa)
      q_sum <- q_sum + sum(meth_a * (total_a - meth_a)) / k_self
      q_den <- q_den + pair_weight
      n_variant_pairs <- n_variant_pairs +
        length(which(wa >= 2)) + length(ia) * (length(ia) - 1) / 2
      add_shared_diagnostic(k_self, pair_weight)
    }

    if (aa == length(groups)) next
    for (bb in (aa + 1L):length(groups)) {
      ib <- groups[[bb]]
      xb <- meth_mat[ib, , drop = FALSE]
      wb <- weights[ib]
      shared <- !is.na(meth_mat[ia[[1L]], ]) &
        !is.na(meth_mat[ib[[1L]], ])
      k_shared <- sum(shared)
      if (k_shared < min_shared_cpg) next

      total_b <- sum(wb)
      pair_weight <- total_a * total_b
      meth_a <- colSums(xa[, shared, drop = FALSE] * wa)
      meth_b <- colSums(xb[, shared, drop = FALSE] * wb)
      q_sum <- q_sum + sum(
        meth_a * (total_b - meth_b) +
          (total_a - meth_a) * meth_b
      ) / k_shared
      q_den <- q_den + pair_weight
      n_variant_pairs <- n_variant_pairs + length(ia) * length(ib)
      add_shared_diagnostic(k_shared, pair_weight)
    }
  }

  shared_cpg_distribution <- if (length(shared_values)) {
    pair_totals <- tapply(shared_weights, shared_values, sum)
    data.frame(
      Shared_CpGs = as.integer(names(pair_totals)),
      Weighted_Eligible_Pairs = as.numeric(pair_totals),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }
  list(
    value = if (q_den > 0) q_sum / q_den else NA_real_,
    status = if (q_den > 0) "ok" else "no_eligible_pairs",
    n_eligible_variant_pairs = as.numeric(n_variant_pairs),
    weighted_eligible_pairs = as.numeric(q_den),
    median_shared_cpg = .panda_weighted_median(shared_values, shared_weights),
    shared_cpg_distribution = shared_cpg_distribution
  )
}

.panda_empty_heterogeneity_scores <- function(status, min_shared_cpg,
                                               min_window_coverage) {
  data.frame(
    Metric = c("Amplicon PDR", "Window Epipolymorphism", "Amplicon qFDRP"),
    Value = rep(NA_real_, 3L),
    Status = rep(status, 3L),
    N_Eligible_Records = c(0, NA, NA),
    Weighted_Eligible_Count = c(0, NA, NA),
    N_Eligible_Windows = c(NA, 0, NA),
    N_Eligible_Variant_Pairs = c(NA, NA, 0),
    Weighted_Eligible_Pairs = c(NA, NA, 0),
    Median_Shared_CpGs = rep(NA_real_, 3L),
    Min_Shared_CpGs = c(NA, NA, min_shared_cpg),
    Min_Window_Coverage = c(NA, min_window_coverage, NA),
    stringsAsFactors = FALSE
  )
}

calculate_heterogeneity <- function(res_obj, min_shared_cpg = 4L,
                                    min_window_coverage = 2L,
                                    cluster_method = "kmeans",
                                    k = 2L,
                                    seed = 11L) {
  min_shared_cpg <- suppressWarnings(as.integer(min_shared_cpg[[1L]]))
  min_window_coverage <- suppressWarnings(as.integer(min_window_coverage[[1L]]))
  if (is.na(min_shared_cpg) || min_shared_cpg < 1L) {
    stop("min_shared_cpg must be a positive integer.", call. = FALSE)
  }
  if (is.na(min_window_coverage) || min_window_coverage < 1L) {
    stop("min_window_coverage must be a positive integer.", call. = FALSE)
  }
  df_long <- res_obj$long_data
  cpg_sites <- res_obj$genome_info$cpg_pos
  empty_scores <- .panda_empty_heterogeneity_scores(
    "no_input_data", min_shared_cpg, min_window_coverage
  )
  if (is.null(df_long) || nrow(df_long) == 0 || length(cpg_sites) == 
      0) 
    return(list(scores = empty_scores, meth_mat = matrix(NA), 
                clusters = NULL, pdr_by_cpg = data.frame(), epipolymorphism_by_window = data.frame(), 
                qfdrp_by_cpg = data.frame(),
                qfdrp_shared_cpg_distribution = data.frame()))
  if (!"Count" %in% names(df_long)) 
    df_long$Count <- 1L
  df_long$Count[is.na(df_long$Count) | df_long$Count < 1] <- 1L
  molecule_counts <- df_long %>% distinct(ReadID, Count) %>% 
    group_by(ReadID) %>% summarise(Count = max(Count), .groups = "drop")
  meth_mat <- df_long %>% select(ReadID, Position, Methylation) %>% 
    distinct(ReadID, Position, .keep_all = TRUE) %>% pivot_wider(names_from = Position, 
                                                                 values_from = Methylation) %>% column_to_rownames("ReadID")
  miss <- setdiff(as.character(cpg_sites), colnames(meth_mat))
  if (length(miss) > 0) 
    for (cc in miss) meth_mat[[cc]] <- NA
  meth_mat <- meth_mat[, as.character(cpg_sites), drop = FALSE]
  molecule_counts <- molecule_counts[match(rownames(meth_mat), 
                                           molecule_counts$ReadID), ]
  molecule_counts$Count[is.na(molecule_counts$Count)] <- 1L
  weights <- molecule_counts$Count
  observed_n <- rowSums(!is.na(meth_mat))
  discordant <- apply(meth_mat, 1, function(x) {
    x <- x[!is.na(x)]
    length(x) >= 4L && length(unique(x)) > 1L
  })
  eligible <- observed_n >= 4L
  pdr_score <- if (any(eligible)) 
    100 * sum(weights[eligible] * discordant[eligible])/sum(weights[eligible])
  else NA_real_
  epi_rows <- list()
  if (length(cpg_sites) >= 4L) {
    for (ww in seq_len(length(cpg_sites) - 3L)) {
      x <- meth_mat[, ww:(ww + 3L), drop = FALSE]
      keep <- complete.cases(x)
      if (sum(weights[keep]) < min_window_coverage || !any(keep)) 
        next
      patterns <- apply(x[keep, , drop = FALSE], 1, paste0, 
                        collapse = "")
      tab <- tapply(weights[keep], patterns, sum)
      pk <- tab/sum(tab)
      epi_rows[[length(epi_rows) + 1L]] <- data.frame(Window = ww, 
                                                      Start_CpG = cpg_sites[ww], End_CpG = cpg_sites[ww + 
                                                                                                       3L], Coverage = sum(weights[keep]), Epipolymorphism = 1 - 
                                                        sum(pk^2))
    }
  }
  epi_df <- if (length(epi_rows)) 
    bind_rows(epi_rows)
  else data.frame()
  epipoly_score <- if (nrow(epi_df)) 
    mean(epi_df$Epipolymorphism)
  else NA_real_
  qfdrp <- .panda_weighted_qfdrp(
    meth_mat,
    weights,
    min_shared_cpg = min_shared_cpg
  )

  # Optional ASM visualization clustering.  This does not alter any of the
  # three heterogeneity metrics; it only supplies a reproducible ordering for
  # the GUI/CLI ASM heatmap.  Rows are clustered on their CpG methylation
  # profiles after column-wise mean imputation for missing calls.
  clusters <- NULL
  cluster_method <- tolower(as.character(cluster_method[[1L]]))
  k <- suppressWarnings(as.integer(k[[1L]]))
  seed <- suppressWarnings(as.integer(seed[[1L]]))
  if (is.na(k) || k < 2L) k <- 2L
  if (is.na(seed)) seed <- 11L
  # stats::kmeans requires more observations than requested centres. When a
  # tiny sample does not satisfy that condition, omit only the optional plot
  # ordering and still return all scientific metrics.
  if (identical(cluster_method, "kmeans") && nrow(meth_mat) > k) {
    cluster_mat <- meth_mat
    for (jj in seq_len(ncol(cluster_mat))) {
      v <- cluster_mat[, jj]
      mu <- mean(v, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0.5
      v[is.na(v)] <- mu
      cluster_mat[, jj] <- v
    }
    if (nrow(unique(cluster_mat)) >= k) {
      had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      km <- stats::kmeans(cluster_mat, centers = k, nstart = 10L)
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
      clusters <- stats::setNames(as.integer(km$cluster), rownames(meth_mat))
    }
  }
  pdr_by_cpg <- lapply(seq_along(cpg_sites), function(cc) {
    keep <- !is.na(meth_mat[, cc]) & eligible
    vals <- meth_mat[keep, , drop = FALSE]
    disc <- if (nrow(vals)) 
      apply(vals, 1, function(x) {
        x <- x[!is.na(x)]
        length(x) >= 4L && length(unique(x)) > 1L
      })
    else logical()
    data.frame(Position = cpg_sites[cc], Coverage = sum(weights[keep]), 
               PDR = if (any(keep)) 
                 sum(weights[keep] * disc)/sum(weights[keep])
               else NA)
  }) %>% bind_rows()
  scores <- data.frame(
    Metric = c("Amplicon PDR", "Window Epipolymorphism", "Amplicon qFDRP"),
    Value = c(pdr_score, epipoly_score, qfdrp$value),
    Status = c(
      if (any(eligible)) "ok" else "no_eligible_records",
      if (nrow(epi_df)) "ok" else "no_eligible_windows",
      qfdrp$status
    ),
    N_Eligible_Records = c(sum(eligible), NA, NA),
    Weighted_Eligible_Count = c(sum(weights[eligible]), NA, NA),
    N_Eligible_Windows = c(NA, nrow(epi_df), NA),
    N_Eligible_Variant_Pairs = c(NA, NA, qfdrp$n_eligible_variant_pairs),
    Weighted_Eligible_Pairs = c(NA, NA, qfdrp$weighted_eligible_pairs),
    Median_Shared_CpGs = c(NA, NA, qfdrp$median_shared_cpg),
    Min_Shared_CpGs = c(NA, NA, min_shared_cpg),
    Min_Window_Coverage = c(NA, min_window_coverage, NA),
    stringsAsFactors = FALSE
  )
  list(scores = scores, meth_mat = meth_mat, clusters = clusters,
       pdr_by_cpg = pdr_by_cpg, epipolymorphism_by_window = epi_df, 
       qfdrp_by_cpg = data.frame(),
       qfdrp_shared_cpg_distribution = qfdrp$shared_cpg_distribution)
}

calculate_quma_stats <- function(res_obj, mode = "Sanger") {
  df_long <- res_obj$long_data
  df_summary <- res_obj$read_summary %>% filter(!str_detect(Pattern, 
                                                            "^excluded"))
  cpg_sites <- res_obj$genome_info$cpg_pos
  if (nrow(df_summary) == 0) 
    return(list(overall = 0, sd_cpg = 0, se_cpg = 0, sd_seq = 0, 
                se_seq = 0, cpg_table = data.frame()))
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
  }
  else {
    overall_meth <- mean(df_long$Methylation) * 100
    cpg_stats <- df_long %>% group_by(Position) %>% summarise(Num_Meth = sum(Methylation), 
                                                              Num_Total = n(), Meth_Pct = mean(Methylation) * 100)
  }
  sd_cpg <- sd(cpg_stats$Meth_Pct, na.rm = T)
  se_cpg <- sd_cpg/sqrt(nrow(cpg_stats))
  sd_seq <- sd(df_summary$Meth_Pct, na.rm = T)
  se_seq <- sd_seq/sqrt(nrow(df_summary))
  full_cpg_stats <- data.frame(Position = cpg_sites) %>% left_join(cpg_stats, 
                                                                   by = "Position") %>% replace_na(list(Num_Meth = 0, Num_Total = 0, 
                                                                                                        Meth_Pct = NA))
  return(list(overall = overall_meth, sd_cpg = sd_cpg, se_cpg = se_cpg, 
              sd_seq = sd_seq, se_seq = se_seq, cpg_table = full_cpg_stats))
}

analyze_group_comparison <- function(res_list_A, res_list_B, 
                                     genome_seq, g1_name = "Group 1", g2_name = "Group 2") {
  if (methods::is(genome_seq, "DNAStringSet")) 
    genome_seq <- genome_seq[[1]]
  cpg_hits <- matchPattern("CG", genome_seq)
  cpg_sites <- start(cpg_hits)
  agg_long <- function(res_list, group_name) {
    do.call(bind_rows, lapply(names(res_list), function(nm) {
      df <- res_list[[nm]]$long_data
      if (is.null(df) || nrow(df) == 0) 
        return(NULL)
      if (!"Count" %in% names(df)) {
        cnt <- str_extract(df$ReadID, "(?<=Count)\\d+")
        df$Count <- ifelse(is.na(cnt), 1, as.integer(cnt))
      }
      df$SampleID <- nm
      df$Group <- group_name
      return(df)
    }))
  }
  df_A <- agg_long(res_list_A, g1_name)
  df_B <- agg_long(res_list_B, g2_name)
  if (is.null(df_A) || is.null(df_B)) 
    return(NULL)
  combined_long <- bind_rows(df_A, df_B)
  agg_site <- function(df) {
    df %>% group_by(Position) %>% summarise(
      Meth = sum(Methylation * Count, na.rm = TRUE),
      Total = sum(Count[is.finite(Methylation)], na.rm = TRUE),
      Pct = ifelse(Total > 0, Meth/Total * 100, NA_real_),
      .groups = "drop"
    )
  }
  site_A <- agg_site(df_A)
  site_B <- agg_site(df_B)
  full_site_table <- data.frame(Position = cpg_sites) %>% left_join(site_A, 
                                                                    by = "Position") %>% rename(Meth_1 = Meth, Total_1 = Total, 
                                                                                                Pct_1 = Pct) %>% left_join(site_B, by = "Position") %>% 
    rename(Meth_2 = Meth, Total_2 = Total, Pct_2 = Pct) %>% 
    replace_na(list(Meth_1 = 0, Total_1 = 0, Meth_2 = 0, Total_2 = 0))
  sample_site <- function(res_list, group_name) {
    bind_rows(lapply(names(res_list), function(nm) {
      rr <- if (!is.null(res_list[[nm]]$metrics)) 
        res_list[[nm]]$metrics
      else res_list[[nm]]
      dd <- rr$long_data
      if (is.null(dd) || !nrow(dd)) 
        return(NULL)
      if (!"Count" %in% names(dd)) 
        dd$Count <- 1L
      dd %>% group_by(Position) %>% summarise(
        Pct = {
          eligible <- is.finite(Methylation)
          if (any(eligible)) {
            sum(Methylation[eligible] * Count[eligible], na.rm = TRUE) /
              sum(Count[eligible], na.rm = TRUE) * 100
          } else NA_real_
        },
        .groups = "drop"
      ) %>%
        mutate(SampleID = nm, Group = group_name)
    }))
  }
  sample_A <- sample_site(res_list_A, g1_name)
  sample_B <- sample_site(res_list_B, g2_name)
  p_vals <- numeric(nrow(full_site_table))
  pooled_p_vals <- numeric(nrow(full_site_table))
  replicate_p_vals <- rep(NA, nrow(full_site_table))
  for (i in 1:nrow(full_site_table)) {
    m1 <- full_site_table$Meth_1[i]
    u1 <- full_site_table$Total_1[i] - m1
    m2 <- full_site_table$Meth_2[i]
    u2 <- full_site_table$Total_2[i] - m2
    if (full_site_table$Total_1[i] == 0 || full_site_table$Total_2[i] ==
        0) {
      pooled_p_vals[i] <- NA
    }
    else {
      mat <- matrix(c(m1, u1, m2, u2), nrow = 2, byrow = TRUE)
      pooled_p_vals[i] <- fisher.test(mat)$p.value
    }
    if (nrow(sample_A) && nrow(sample_B)) {
      a <- sample_A$Pct[sample_A$Position == cpg_sites[i]]
      b <- sample_B$Pct[sample_B$Position == cpg_sites[i]]
      if (sum(is.finite(a)) >= 2L && sum(is.finite(b)) >= 
          2L) 
        replicate_p_vals[i] <- tryCatch(t.test(a, b)$p.value, 
                                        error = function(e) NA)
    }
  }
  p_vals <- ifelse(is.finite(replicate_p_vals), replicate_p_vals, 
                   pooled_p_vals)
  full_site_table$P_Value <- p_vals
  full_site_table$FDR <- p.adjust(p_vals, method = "BH")
  full_site_table$P_Value_Source <- ifelse(
    is.finite(replicate_p_vals), "sample_level_welch_t",
    ifelse(is.finite(pooled_p_vals), "pooled_read_level_fisher", "not_estimable")
  )
  # Do not substitute pooled reads or dereplicated variants for biological
  # samples in the overall group test. Such a test would constitute
  # pseudoreplication and is deliberately not performed.
  pooled_u_test_p <- NA_real_
  sample_overall <- function(x) {
    vapply(x, function(r) {
      rr <- if (!is.null(r$metrics)) 
        r$metrics
      else r
      d <- rr$long_data
      if (!nrow(d)) 
        return(NA)
      eligible <- is.finite(d$Methylation)
      if (!any(eligible)) return(NA_real_)
      if (!"Count" %in% names(d)) d$Count <- 1
      sum(d$Methylation[eligible] * d$Count[eligible], na.rm = TRUE) /
        sum(d$Count[eligible], na.rm = TRUE) * 100
    }, numeric(1))
  }
  oa <- sample_overall(res_list_A)
  ob <- sample_overall(res_list_B)
  oa <- oa[is.finite(oa)]
  ob <- ob[is.finite(ob)]
  sample_level_test_available <- length(oa) >= 2L && length(ob) >= 2L
  u_test_p <- if (sample_level_test_available)
    tryCatch(wilcox.test(oa, ob, exact = FALSE)$p.value, 
             error = function(e) NA)
  else NA_real_
  sample_values <- data.frame(
    Sample = c(names(oa), names(ob)),
    Group = c(rep(g1_name, length(oa)), rep(g2_name, length(ob))),
    Overall_Methylation = c(unname(oa), unname(ob)),
    stringsAsFactors = FALSE
  )
  sum_df <- data.frame(
    Group = c(g1_name, g2_name),
    N = c(length(oa), length(ob)),
    Mean = c(mean(oa), mean(ob)),
    SD = c(sd(oa), sd(ob)),
    stringsAsFactors = FALSE
  )
  sum_df$Mean[!is.finite(sum_df$Mean)] <- NA_real_
  sum_df$SD[!is.finite(sum_df$SD)] <- NA_real_
  return(list(site_table = full_site_table, combined_long = combined_long, 
              u_test_p = u_test_p, pooled_u_test_p = pooled_u_test_p, 
              replicate_p_used = any(is.finite(replicate_p_vals)), 
              overall_test = if (sample_level_test_available)
                "sample_level_wilcoxon" else "not_estimable",
              pooled_overall_test = "not_performed_to_avoid_pseudoreplication",
              site_test_uses_replicates = any(is.finite(replicate_p_vals)),
              sample_values = sample_values,
              sample_summary = sum_df,
              summary = sum_df, g1_name = g1_name,
              g2_name = g2_name))
}
