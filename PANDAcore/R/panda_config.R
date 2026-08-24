# PANDA configuration helpers

#' Read a PANDA YAML or JSON configuration file
#'
#' @param path Path to a `.yml`, `.yaml`, or `.json` file.
#' @return A named list.
#' @export
panda_read_config <- function(path) {
  if (length(path) != 1L || !is.character(path) || !file.exists(path)) {
    stop("Configuration file does not exist: ", path, call. = FALSE)
  }
  
  ext <- tolower(tools::file_ext(path))
  
  cfg <- switch(
    ext,
    
    yml = {
      if (!requireNamespace("yaml", quietly = TRUE)) {
        stop(
          "Package 'yaml' is required for .yml configuration files.",
          call. = FALSE
        )
      }
      yaml::read_yaml(path)
    },
    
    yaml = {
      if (!requireNamespace("yaml", quietly = TRUE)) {
        stop(
          "Package 'yaml' is required for .yaml configuration files.",
          call. = FALSE
        )
      }
      yaml::read_yaml(path)
    },
    
    json = {
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop(
          "Package 'jsonlite' is required for .json configuration files.",
          call. = FALSE
        )
      }
      jsonlite::fromJSON(path, simplifyVector = FALSE)
    },
    
    stop(
      "Configuration file must have a .yml, .yaml, or .json extension.",
      call. = FALSE
    )
  )
  
  if (!is.list(cfg) || is.null(names(cfg))) {
    stop(
      "The configuration file must contain a named mapping/object.",
      call. = FALSE
    )
  }
  
  cfg
}


#' Validate a PANDA configuration
#'
#' @param config A named list.
#' @return The validated configuration, invisibly.
#' @export
panda_validate_config <- function(config) {
  if (!is.list(config) || is.null(names(config))) {
    stop("config must be a named list.", call. = FALSE)
  }
  
  required <- c("mode", "input", "output")
  missing <- setdiff(required, names(config))
  
  if (length(missing) > 0L) {
    stop(
      "Missing required configuration field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  mode <- tolower(as.character(config$mode[[1L]]))
  
  if (!mode %in% c("sanger", "amplicon", "ngs")) {
    stop(
      "mode must be one of: sanger, amplicon, ngs.",
      call. = FALSE
    )
  }
  
  input_path <- as.character(config$input[[1L]])
  
  if (!file.exists(input_path) && !dir.exists(input_path)) {
    stop("Input path does not exist: ", input_path, call. = FALSE)
  }
  
  output_path <- as.character(config$output[[1L]])
  
  if (length(output_path) != 1L || !nzchar(output_path)) {
    stop("output must be a non-empty path.", call. = FALSE)
  }
  
  invisible(config)
}


#' List supported PANDA sequence input files
#'
#' @param path A file or directory.
#' @return A character vector of input paths.
#' @export
panda_input_files <- function(path) {
  if (length(path) != 1L || !is.character(path)) {
    stop("path must be a single character string.", call. = FALSE)
  }
  
  # A directory may also return TRUE for file.exists() on some systems.
  if (file.exists(path) && !dir.exists(path)) {
    return(normalizePath(path, mustWork = TRUE))
  }
  
  if (!dir.exists(path)) {
    stop("Input path does not exist: ", path, call. = FALSE)
  }
  
  files <- list.files(
    path = path,
    pattern = "\\.(ab1|abi|fasta?|fastq)(\\.gz)?$",
    ignore.case = TRUE,
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(files) == 0L) {
    stop("No supported input files were found in: ", path, call. = FALSE)
  }
  
  sort(normalizePath(files, mustWork = TRUE))
}