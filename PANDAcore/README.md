# PANDAcore

`PANDAcore` is the shared R package layer for the PANDA platform. It contains
the non-UI functions used by the command-line workflow, including:

- bisulfite-aware global alignment through `pwalign`;
- FASTA/FASTQ/AB1 input handling and read dereplication;
- methylation statistics and CpG-level summaries; and
- PDR, window epipolymorphism, qFDRP, and ASM clustering calculations.

The heterogeneity implementation treats Sanger records as equally weighted
clones and NGS `Count` values as observed retained-read abundance. Amplicon
qFDRP is exactly equivalent to expanding dereplicated variants by `Count`, so
same-variant read pairs are retained with zero distance. Ineligible metrics are
reported as `NA` together with eligibility diagnostics and thresholds.

The Shiny application and CLI are the user-facing interfaces. The package is
kept in the repository so that the computational layer can be tested and
versioned independently from the GUI.

For local development, install it into the project environment from the
repository root:

```r
renv::install("./PANDAcore", rebuild = TRUE)
```

The public entry points are listed in `PANDAcore/NAMESPACE`. The CLI and the
local Shiny application load this package through the project environment;
install it once from the repository root before running either interface.
