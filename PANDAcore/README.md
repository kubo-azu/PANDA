# PANDAcore

This is the shared R package layer for PANDA. The Shiny application and the
command-line interface will call the same functions from this package so that
identical inputs and parameters produce identical analysis results.

The package is under active development. The first milestone is to move the
tested, non-UI functions from `app.R` into `R/` without changing their behavior.
