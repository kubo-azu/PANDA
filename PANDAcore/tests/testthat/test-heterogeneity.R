make_heterogeneity_result <- function(mat, weights = rep.int(1L, nrow(mat))) {
  positions <- seq_len(ncol(mat)) * 10L
  rows <- lapply(seq_len(nrow(mat)), function(ii) {
    keep <- !is.na(mat[ii, ])
    data.frame(
      ReadID = rownames(mat)[[ii]],
      Position = positions[keep],
      Methylation = as.numeric(mat[ii, keep]),
      Count = weights[[ii]],
      stringsAsFactors = FALSE
    )
  })
  list(
    long_data = dplyr::bind_rows(rows),
    genome_info = list(cpg_pos = positions)
  )
}

test_that("qFDRP includes zero-distance pairs from the same variant", {
  mat <- rbind(A = c(0, 0, 0, 0), B = c(1, 1, 1, 1))
  result <- calculate_heterogeneity(
    make_heterogeneity_result(mat, c(3L, 1L))
  )
  q <- result$scores[result$scores$Metric == "Amplicon qFDRP", ]
  expect_equal(q$Value, 0.5)
  expect_equal(q$Weighted_Eligible_Pairs, 6)
  expect_identical(q$Status, "ok")
  expect_equal(
    sum(result$qfdrp_shared_cpg_distribution$Weighted_Eligible_Pairs), 6
  )
})

test_that("zero heterogeneity is distinct from an inestimable value", {
  mat <- rbind(A = c(0, 0, 0, 0))
  replicated <- calculate_heterogeneity(
    make_heterogeneity_result(mat, 3L)
  )
  lone <- calculate_heterogeneity(
    make_heterogeneity_result(mat, 1L)
  )
  expect_equal(replicated$scores$Value[[3L]], 0)
  expect_identical(replicated$scores$Status[[3L]], "ok")
  expect_true(is.na(lone$scores$Value[[3L]]))
  expect_identical(lone$scores$Status[[3L]], "no_eligible_pairs")
})

test_that("eligibility failures and empty input return NA", {
  partial <- rbind(A = c(0, 0, NA, NA), B = c(1, 1, NA, NA))
  result <- calculate_heterogeneity(
    make_heterogeneity_result(partial), min_shared_cpg = 4L
  )
  expect_true(is.na(result$scores$Value[[3L]]))
  expect_identical(result$scores$Status[[3L]], "no_eligible_pairs")

  empty <- calculate_heterogeneity(list(
    long_data = data.frame(), genome_info = list(cpg_pos = integer())
  ))
  expect_true(all(is.na(empty$scores$Value)))
  expect_true(all(empty$scores$Status == "no_input_data"))
})

make_group_sample <- function(id, methylation, positions = c(1L, 3L)) {
  long_data <- data.frame(
    ReadID = rep(id, length(positions)),
    Position = positions,
    Methylation = methylation,
    Count = 1L,
    stringsAsFactors = FALSE
  )
  list(
    long_data = long_data,
    read_summary = data.frame(
      ReadID = id,
      Meth_Pct = mean(methylation) * 100,
      Pattern = "Passed",
      Count = 1L,
      stringsAsFactors = FALSE
    )
  )
}

test_that("group summaries use the same call-level overall methylation definition", {
  group_a <- list(
    A1 = make_group_sample("A1", c(0, 0)),
    A2 = make_group_sample("A2", c(0, 1)),
    A3 = make_group_sample("A3", c(1, 0))
  )
  group_b <- list(
    B1 = make_group_sample("B1", c(1, 1)),
    B2 = make_group_sample("B2", c(1, 0)),
    B3 = make_group_sample("B3", c(0, 1))
  )
  comparison <- analyze_group_comparison(
    group_a, group_b, Biostrings::DNAString("CGCG"), "A", "B"
  )

  expect_equal(
    comparison$sample_values$Overall_Methylation,
    c(0, 50, 50, 100, 50, 50)
  )
  expect_equal(comparison$sample_summary$Mean, c(100 / 3, 200 / 3))
  expect_equal(comparison$sample_summary$SD, c(50 / sqrt(3), 50 / sqrt(3)))
  expect_identical(
    comparison$overall_test,
    "sample_level_wilcoxon_exact"
  )
  expect_true(comparison$overall_test_exact)
  expect_identical(
    unname(comparison$overall_test_r_method),
    "Wilcoxon rank sum exact test"
  )
  expect_true(is.finite(comparison$u_test_p))
  expect_true(all(comparison$site_table$P_Value_Source == "sample_level_welch_t"))
})

test_that("small comparisons without ties use the R exact p-value", {
  positions <- seq.int(1L, 19L, by = 2L)
  make_proportion_sample <- function(id, methylated) {
    make_group_sample(
      id,
      c(rep.int(1, methylated), rep.int(0, 10L - methylated)),
      positions
    )
  }
  group_a <- list(
    A1 = make_proportion_sample("A1", 0L),
    A2 = make_proportion_sample("A2", 1L),
    A3 = make_proportion_sample("A3", 2L)
  )
  group_b <- list(
    B1 = make_proportion_sample("B1", 8L),
    B2 = make_proportion_sample("B2", 9L),
    B3 = make_proportion_sample("B3", 10L)
  )
  comparison <- analyze_group_comparison(
    group_a, group_b,
    Biostrings::DNAString(paste(rep("CG", 10L), collapse = "")),
    "A", "B"
  )
  expect_equal(comparison$u_test_p, 0.1)
  expect_identical(comparison$overall_test, "sample_level_wilcoxon_exact")
  expect_true(comparison$overall_test_exact)
  expect_match(comparison$overall_test_r_method, "exact")
})

test_that("R selects the asymptotic approximation for larger comparisons", {
  comparison <- .panda_sample_level_wilcoxon(1:50, 51:100)
  expect_identical(comparison$method, "sample_level_wilcoxon_asymptotic")
  expect_false(comparison$exact)
  expect_identical(
    unname(comparison$r_method),
    "Wilcoxon rank sum test with continuity correction"
  )
  expect_true(is.finite(comparison$p.value))
})

test_that("insufficient biological replication is not replaced by read-level inference", {
  deep_sample <- function(id, values) {
    x <- make_group_sample(id, values)
    x$long_data$Count <- 10000L
    x
  }
  comparison <- analyze_group_comparison(
    list(A1 = deep_sample("A1", c(0, 0))),
    list(B1 = deep_sample("B1", c(1, 1))),
    Biostrings::DNAString("CGCG"), "A", "B"
  )
  expect_true(is.na(comparison$u_test_p))
  expect_identical(comparison$overall_test, "not_estimable")
  expect_true(is.na(comparison$pooled_u_test_p))
  expect_identical(
    comparison$pooled_overall_test,
    "not_performed_to_avoid_pseudoreplication"
  )
  expect_true(all(
    comparison$site_table$P_Value_Source == "not_estimable"
  ))
  expect_true(all(is.na(comparison$site_table$P_Value)))
  expect_true(all(is.na(comparison$site_table$FDR)))
  expect_true(all(comparison$site_table$Test_Status == "insufficient_samples"))
  expect_equal(comparison$site_table$N_Samples_1, c(1L, 1L))
  expect_equal(comparison$site_table$N_Samples_2, c(1L, 1L))
  expect_equal(comparison$site_table$Pct_1, c(0, 0))
  expect_equal(comparison$site_table$Pct_2, c(100, 100))

})

test_that("zero coverage is represented as missing rather than unmethylated", {
  comparison <- analyze_group_comparison(
    list(A1 = make_group_sample("A1", c(0, 1))),
    list(B1 = make_group_sample("B1", 1, positions = 1L)),
    Biostrings::DNAString("CGCG"), "A", "B"
  )
  missing_row <- comparison$site_table[
    comparison$site_table$Position == 3L, , drop = FALSE
  ]
  expect_true(is.na(missing_row$Pct_2))
  expect_true(is.na(missing_row$P_Value))
  expect_identical(missing_row$P_Value_Source, "not_estimable")
})


test_that("Welch failures remain missing and BH uses only estimable sites", {
  positions <- c(1L, 3L, 5L, 7L)
  make_group <- function(prefix, values) {
    setNames(lapply(seq_len(nrow(values)), function(i)
      make_group_sample(paste0(prefix, i), values[i, ], positions)),
      paste0(prefix, seq_len(nrow(values))))
  }
  a <- rbind(c(0, 0, 0, 0), c(1, 0, 0, 1), c(0, 0, 1, 0))
  b <- rbind(c(1, 1, 1, NA), c(1, 1, 1, NA), c(0, 1, 0, NA))
  result <- analyze_group_comparison(make_group("A", a), make_group("B", b),
                                     Biostrings::DNAString("CGCGCGCG"))
  site <- result$site_table
  expected <- c(t.test(a[, 1] * 100, b[, 1] * 100)$p.value,
                t.test(a[, 3] * 100, b[, 3] * 100)$p.value)
  expect_equal(site$P_Value[c(1, 3)], expected)
  expect_equal(site$FDR[c(1, 3)], p.adjust(expected, method = "BH"))
  expect_true(all(is.na(site$P_Value[c(2, 4)])))
  expect_true(all(is.na(site$FDR[c(2, 4)])))
  expect_equal(site$Test_Status, c("ok", "zero_variance_in_both_groups", "ok",
                                  "no_coverage_in_one_or_both_groups"))
  expect_equal(site$N_Samples_2, c(3L, 3L, 3L, 0L))
})

test_that("one constant group still permits an estimable Welch test", {
  result <- analyze_group_comparison(
    list(A1 = make_group_sample("A1", c(0, 0)),
         A2 = make_group_sample("A2", c(0, 0))),
    list(B1 = make_group_sample("B1", c(0, 0)),
         B2 = make_group_sample("B2", c(1, 1))),
    Biostrings::DNAString("CGCG"))
  expect_true(all(is.finite(result$site_table$P_Value)))
  expect_true(all(result$site_table$Test_Status == "ok"))
})

test_that("qFDRP does not distinguish phasing with identical site frequencies", {
  bimodal <- rbind(A = c(0, 0, 0, 0), B = c(1, 1, 1, 1))
  diverse <- as.matrix(expand.grid(rep(list(0:1), 4)))
  rownames(diverse) <- paste0("V", seq_len(nrow(diverse)))
  x <- calculate_heterogeneity(make_heterogeneity_result(bimodal, c(80L, 80L)))
  y <- calculate_heterogeneity(make_heterogeneity_result(diverse, rep(10L, 16)))
  expect_equal(x$scores$Value[3], 160 / 159 * 0.5)
  expect_equal(x$scores$Value[3], y$scores$Value[3])
  expect_gt(y$scores$Value[1], x$scores$Value[1])
  expect_gt(y$scores$Value[2], x$scores$Value[2])
})

test_that("shared CpG minimum filters pairs without truncating their distance", {
  mat <- rbind(A = c(0, 0, 0, 0), B = c(0, 0, 1, 1), C = c(1, 1, NA, NA))
  input <- make_heterogeneity_result(mat)
  low <- calculate_heterogeneity(input, min_shared_cpg = 2L)
  high <- calculate_heterogeneity(input, min_shared_cpg = 4L)
  none <- calculate_heterogeneity(input, min_shared_cpg = 5L)
  expect_equal(low$scores$Value[3], (0.5 + 1 + 1) / 3)
  expect_equal(high$scores$Value[3], 0.5)
  expect_equal(low$scores$Weighted_Eligible_Pairs[3], 3)
  expect_equal(high$scores$Weighted_Eligible_Pairs[3], 1)
  expect_true(is.na(none$scores$Value[3]))
  full <- make_heterogeneity_result(mat[1:2, ])
  expect_equal(calculate_heterogeneity(full, min_shared_cpg = 2L)$scores$Value[3],
               calculate_heterogeneity(full, min_shared_cpg = 4L)$scores$Value[3])
})
