# Configuration helpers shared by the GUI and CLI.

#' Read a PANDA YAML configuration file.
#'
#' @param path Path to a YAML file.
#' @return A named list.
#' @export
panda_read_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read configuration files.", call. = FALSE)
  }
  if (length(path) != 1L || !file.exists(path)) {
    stop("Configuration file does not exist: ", path, call. = FALSE)
  }
  cfg <- yaml::read_yaml(path)
  if (!is.list(cfg)) {
    stop("The configuration file must contain a YAML mapping.", call. = FALSE)
  }
  cfg
}

#' Validate the common PANDA configuration fields.
#'
#' @param config A named list.
#' @return The validated configuration, invisibly.
#' @export
panda_validate_config <- function(config) {
  if (!is.list(config)) {
    stop("config must be a named list.", call. = FALSE)
  }
  required <- c("mode", "input", "output")
  missing <- setdiff(required, names(config))
  if (length(missing)) {
    stop("Missing required configuration field(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  mode <- tolower(as.character(config$mode[[1L]]))
  if (!mode %in% c("sanger", "amplicon", "ngs")) {
    stop("mode must be one of: sanger, amplicon, ngs.", call. = FALSE)
  }
  if (!dir.exists(config$input) && !file.exists(config$input)) {
    stop("Input path does not exist: ", config$input, call. = FALSE)
  }
  if (length(config$output) != 1L || !nzchar(config$output)) {
    stop("output must be a non-empty path.", call. = FALSE)
  }
  invisible(config)
}

#' List supported sequence input files.
#'
#' @param path A file or directory.
#' @return A character vector of input paths.
#' @export
panda_input_files <- function(path) {
  if (length(path) != 1L) stop("path must have length one.", call. = FALSE)
  if (file.exists(path)) return(normalizePath(path, mustWork = TRUE))
  if (!dir.exists(path)) stop("Input path does not exist: ", path, call. = FALSE)
  files <- list.files(
    path,
    pattern = "\\.(ab1|abi|fasta?|fastq)(\\.gz)?$",
    ignore.case = TRUE,
    full.names = TRUE,
    recursive = TRUE
  )
  sort(normalizePath(files, mustWork = TRUE))
}
