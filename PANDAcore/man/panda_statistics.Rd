\name{calculate_heterogeneity}
\alias{calculate_heterogeneity}
\alias{calculate_quma_stats}
\alias{analyze_group_comparison}
\title{Summarize Methylation and Compare Groups}
\usage{
calculate_heterogeneity(
  res_obj,
  min_shared_cpg = 4L,
  min_window_coverage = 2L,
  cluster_method = "kmeans",
  k = 2L,
  seed = 11L
)

calculate_quma_stats(res_obj, mode = "Sanger")

analyze_group_comparison(
  res_list_A,
  res_list_B,
  genome_seq,
  g1_name = "Group 1",
  g2_name = "Group 2"
)
}
\arguments{
  \item{res_obj}{A PANDA alignment result.}
  \item{min_shared_cpg}{Minimum shared CpGs required for an eligible qFDRP pair.}
  \item{min_window_coverage}{Minimum read-count-weighted coverage for a four-CpG window.}
  \item{cluster_method}{Optional display-clustering method; currently \code{"kmeans"}.}
  \item{k}{Number of display clusters.}
  \item{seed}{Random seed used for display clustering.}
  \item{mode}{Either \code{"Sanger"} or \code{"NGS"}; NGS summaries use retained read-count weights.}
  \item{res_list_A, res_list_B}{Named lists of independent sample results for two groups.}
  \item{genome_seq}{Reference sequence as a \code{DNAString} or one-record \code{DNAStringSet}.}
  \item{g1_name, g2_name}{Labels for the two groups.}
}
\value{
\code{calculate_heterogeneity()} returns metric scores and diagnostic tables.
\code{calculate_quma_stats()} returns overall and CpG-level methylation
summaries. \code{analyze_group_comparison()} returns sample-level overall
comparisons and replicate-level CpG comparisons with test-status fields.
}
\description{
Computes PANDA methylation summaries, amplicon-adapted heterogeneity metrics,
and two-group comparisons while retaining observed NGS read-count weights.
}
