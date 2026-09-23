\name{panda_dereplicate_reads}
\alias{panda_dereplicate_reads}
\title{Dereplicate Identical Amplicon Reads}
\usage{
panda_dereplicate_reads(reads, min_count = 1L, max_unique_reads = Inf)
}
\arguments{
  \item{reads}{A \code{DNAStringSet} containing input reads.}
  \item{min_count}{Minimum observed multiplicity required to retain a sequence variant.}
  \item{max_unique_reads}{Maximum number of sequence variants to retain, or \code{Inf}.}
}
\value{
A list containing representative reads, their integer counts, the input read
count, and the number of retained unique sequences.
}
\description{
Collapses reads that are identical across the complete retained sequence.
Representatives are ordered by decreasing observed read count.
}
