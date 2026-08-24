# PANDAcore

`PANDAcore` is the shared R package layer for the PANDA platform. It contains
the non-UI functions used by the command-line workflow, including:

- bisulfite-aware global alignment through `pwalign`;
- FASTA/FASTQ/AB1 input handling and read dereplication;
- methylation statistics and CpG-level summaries; and
- PDR, window epipolymorphism, qFDRP, and ASM clustering calculations.

The Shiny application and CLI are the user-facing interfaces. The package is
kept in the repository so that the computational layer can be tested and
versioned independently from the GUI.

For local development, install it into the project environment from the
repository root:

```r
renv::install("./PANDAcore", rebuild = TRUE)
```

The public entry points are listed in `PANDAcore/NAMESPACE`. The CLI scripts
load this package through the project environment; end users normally do not
need to install `PANDAcore` separately.
