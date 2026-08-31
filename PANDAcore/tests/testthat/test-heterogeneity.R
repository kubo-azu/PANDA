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
  expect_identical(comparison$overall_test, "sample_level_wilcoxon")
  expect_true(is.finite(comparison$u_test_p))
  expect_true(all(comparison$site_table$P_Value_Source == "sample_level_welch_t"))
})

test_that("insufficient biological replication is not replaced by read-level inference", {
  comparison <- analyze_group_comparison(
    list(A1 = make_group_sample("A1", c(0, 0))),
    list(B1 = make_group_sample("B1", c(1, 1))),
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
    comparison$site_table$P_Value_Source == "pooled_read_level_fisher"
  ))
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
