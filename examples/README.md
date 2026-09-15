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

## Scope of these examples

Experiments A-C demonstrate PANDA features using synthetic data. The example
configurations support analysis and plotting of the bundled data; they do not
constitute a general robustness benchmark. Manuscript-specific figure and
summary scripts are not included in this public repository.

If more than one R installation is present, set `R_PANDA_RSCRIPT` to the
intended `Rscript` executable when using the installed `panda` launcher.
