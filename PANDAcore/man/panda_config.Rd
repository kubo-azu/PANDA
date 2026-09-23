\name{panda_read_config}
\alias{panda_read_config}
\alias{panda_validate_config}
\alias{panda_input_files}
\title{Read and Validate PANDA Inputs}
\usage{
panda_read_config(path)

panda_validate_config(config)

panda_input_files(path)
}
\arguments{
  \item{path}{Path to a PANDA JSON/YAML configuration, a supported sequence
  file, or a directory containing supported sequence files.}
  \item{config}{A named list containing PANDA configuration values.}
}
\value{
\code{panda_read_config()} returns a named list. \code{panda_validate_config()}
returns the validated configuration invisibly. \code{panda_input_files()}
returns normalized paths to supported input files.
}
\description{
Utilities used by the command-line workflow to read configuration files,
validate required settings, and discover AB1, FASTA, and FASTQ inputs.
}
