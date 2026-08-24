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
