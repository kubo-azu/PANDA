# PANDA example configurations

These files are small, reproducible configuration examples for the command-line
workflow. They are not analysis results.

Run an example from the repository root after restoring the `renv` environment:

```bash
panda analyze examples/analysis_C.json
panda plot results/demo_C_group
```

The examples use the synthetic data under `PANDA_Demo/`:

- `analysis_C.json`: NGS WT/KO group comparison.
- `analysis_B.json`: mixed-haplotype NGS data without motif filtering.
- `analysis_B_motif_A.json` and `analysis_B_motif_B.json`: motif-filtered
  haplotype analyses.
- `analysis_sanger_A.json`, `analysis_sanger_B.json`, and
  `analysis_sanger_C.json`: Sanger examples.

Motif files under `examples/motifs/` contain one literal DNA motif per line.
Paths in these configurations are relative to the repository root.

Group membership files under `examples/groups/` contain one sample ID per
line. They are used by `panda compare` and avoid embedding long sample lists in
the command itself.

## Manuscript demo re-analysis

The complete, option-explicit manuscript workflow is available as:

```bash
bash tools/reanalyze_manuscript_demos.sh
```

It analyzes all reads in Experiments A-C with both NGS and Sanger modes,
performs the Experiment B motif-filtered runs, creates plots and group
comparisons, summarizes the synthetic ground-truth agreement, and runs the
qFDRP shared-CpG sensitivity analysis. Results are written to
`results/manuscript_reanalysis`.

Ground-truth tables report methylation error and the retained-record recovery
rate. Accuracy and qFDRP sensitivity PDFs are written separately for NGS and
Sanger so the higher-depth NGS validation can be presented independently while
retaining the Sanger validation as a complementary result.

The script refuses to mix a new run with an existing non-empty result
directory. A different destination or worker count can be specified without
editing the script:

```bash
PANDA_MANUSCRIPT_OUTPUT=results/manuscript_reanalysis_v2 \
PANDA_MANUSCRIPT_WORKERS=16 \
bash tools/reanalyze_manuscript_demos.sh
```

If more than one R installation is present, set `R_PANDA_RSCRIPT` to the
intended `Rscript` executable. The manuscript script and the installed
`panda` launcher will then use the same R runtime.
