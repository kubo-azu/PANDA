#' Collapse identical reads before alignment.
#'
#' Reads are sorted by descending multiplicity. The returned representative
#' names are stable and encode the abundance in a separate Count vector.
#'
#' @param reads A Biostrings DNAStringSet.
#' @param min_count Minimum abundance required for a representative.
#' @param max_unique_reads Optional cap on the number of representatives.
#' @return A list containing representative reads, counts, and read totals.
#' @export
panda_dereplicate_reads <- function(reads, min_count = 1L,
                                    max_unique_reads = Inf) {
  if (!inherits(reads, "DNAStringSet")) {
    stop("reads must be a Biostrings DNAStringSet.", call. = FALSE)
  }
  
  if (length(reads) == 0L) {
    return(list(
      reads = Biostrings::DNAStringSet(),
      counts = integer(),
      raw_n = 0L,
      unique_n = 0L
    ))
  }
  
  sequence_strings <- as.character(reads)
  counts <- sort(table(sequence_strings), decreasing = TRUE)
  
  min_count <- as.integer(min_count)
  if (is.na(min_count) || min_count < 1L) {
    stop("min_count must be a positive integer.", call. = FALSE)
  }
  counts <- counts[counts >= min_count]
  
  if (is.finite(max_unique_reads)) {
    max_unique_reads <- as.integer(max_unique_reads)
    if (max_unique_reads < 1L) {
      stop("max_unique_reads must be positive or Inf.", call. = FALSE)
    }
    counts <- utils::head(counts, max_unique_reads)
  }
  
  representative_reads <- Biostrings::DNAStringSet(names(counts))
  names(representative_reads) <- paste0(
    "Rank", seq_along(counts),
    "_Count", as.integer(counts)
  )
  
  list(
    reads = representative_reads,
    counts = as.integer(counts),
    raw_n = length(sequence_strings),
    unique_n = length(representative_reads)
  )
}
