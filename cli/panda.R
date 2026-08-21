#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(paste0(
    "PANDA command-line interface (development version)\n\n",
    "Usage:\n",
    "  Rscript cli/panda.R --config analysis.yml\n",
    "  Rscript cli/panda.R --help\n"
  ))
}

if (!length(args) || "--help" %in% args || "-h" %in% args) {
  usage()
  quit(status = 0L)
}

if (!"--config" %in% args) {
  usage()
  stop("--config is required in this development release.", call. = FALSE)
}

idx <- match("--config", args)
if (idx == length(args)) stop("--config requires a path.", call. = FALSE)
config_path <- args[[idx + 1L]]

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "cli/panda.R"
pkg_path <- file.path(dirname(normalizePath(script_path, mustWork = TRUE)), "..", "PANDAcore")
if (!dir.exists(pkg_path)) stop("PANDAcore directory was not found.", call. = FALSE)

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Install the 'yaml' package before using --config.", call. = FALSE)
}

source(file.path(pkg_path, "R", "panda_config.R"), local = TRUE)
cfg <- panda_read_config(config_path)
panda_validate_config(cfg)
cat("Configuration validated successfully. Analysis engine wiring is the next milestone.\n")
