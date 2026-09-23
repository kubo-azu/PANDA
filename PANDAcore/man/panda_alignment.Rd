\name{run_bisulfite_alignment}
\alias{run_bisulfite_alignment}
\alias{process_ab1_files}
\title{Prepare and Align Bisulfite Amplicon Sequences}
\usage{
run_bisulfite_alignment(
  genome_seq,
  reads_set,
  min_identity = 90,
  min_conversion = 95,
  return_alignments = TRUE,
  workers = 16L
)

process_ab1_files(
  file_paths,
  file_names,
  trim_start = 20,
  trim_end = 20
)
}
\arguments{
  \item{genome_seq}{Reference amplicon sequence as a \code{DNAString} or a one-record \code{DNAStringSet}.}
  \item{reads_set}{A named \code{DNAStringSet} containing reads or dereplicated sequence variants.}
  \item{min_identity}{Minimum alignment identity percentage.}
  \item{min_conversion}{Minimum non-CpG conversion percentage.}
  \item{return_alignments}{Whether to retain printable pairwise alignments.}
  \item{workers}{Number of parallel workers, capped at 16.}
  \item{file_paths}{Character vector of AB1 file paths.}
  \item{file_names}{Names assigned to the resulting sequences.}
  \item{trim_start}{Number of bases trimmed from the start of each AB1-derived sequence.}
  \item{trim_end}{Number of bases trimmed from the end of each AB1-derived sequence.}
}
\value{
\code{run_bisulfite_alignment()} returns reference-coordinate CpG calls,
per-record quality summaries, reference information, counts, and optional
alignment text. \code{process_ab1_files()} returns a \code{DNAStringSet}.
}
\description{
Aligns bisulfite-converted amplicon sequences to a defined reference and calls
methylation at CpG positions. AB1 chromatograms can be converted to trimmed
primary sequences before alignment.
}
