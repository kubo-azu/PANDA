# 🐼 PANDA: Phased ANalysis of DNA Amplicons 🧬

**Current release: 1.0.0**

**PANDA (Phased ANalysis of DNA Amplicons)** is an R-based platform for
phased DNA-methylation analysis from targeted bisulfite sequencing data. It
combines a Shiny graphical interface with a reproducible command-line
interface, using the same analysis concepts for Sanger and amplicon-NGS data.

**💡 Web version:** PANDA can be explored in a browser at
<https://huggingface.co/spaces/kubo-azu/PANDA>.

**💡 Please cite:** Kubota, A., Kobayashi, H., & Tajima, A. (2026). PANDA:
Read-Level Phased Analysis of DNA Amplicons for Methylation Studies.
*bioRxiv*. <https://doi.org/10.64898/2026.04.01.715790>.

PANDA is designed for questions in which the methylation state of individual
DNA molecules matters, including allele-specific methylation (ASM), imprinting
and loss of imprinting (LOI), haplotype-resolved methylation, and within-sample
heterogeneity.

## 🌟 What PANDA provides

- **Sanger mode** for FASTA and AB1-derived clonal sequences.
- **Amplicon-NGS mode** for FASTA/FASTQ data, dereplication, read-count
  weighting, and optional paired-end handling.
- **Global bisulfite-aware alignment** using `pwalign`, with explicit identity
  and conversion filters.
- **Phased methylation views** including lollipop plots, ASM profiles, and a
  count-weighted binned abundance heatmap.
- **Heterogeneity metrics** including PDR, window epipolymorphism, and qFDRP.
- **In-silico haplotype filtering** from motifs supplied directly or through a
  one-motif-per-line text file.
- **Batch summaries and group comparison** with per-CpG differences and
  multiple-testing correction.
- **GUI and CLI workflows** that write machine-readable tables and PDFs.

PANDA is intended for targeted amplicons, not whole-genome alignment. Read
quality trimming, adapter/primer removal, and other experiment-specific QC
should be performed before PANDA when appropriate.

## 🧩 Platform layout

```text
app.R                 Current Shiny application
cli/panda.R           analysis command
cli/panda_plot.R      plotting command
cli/panda_compare.R   two-group comparison command
PANDAcore/            shared R package layer
bin/panda             user-facing command launcher
install_panda.sh      user-level launcher installer
examples/             reproducible configuration examples
PANDA_Demo/           synthetic demonstration data
```

The CLI writes an analysis bundle containing `PANDA_analysis.json`, summary
tables, per-sample statistics, alignments, and a `plots/` directory. These are
analysis outputs and should normally be kept outside version control.

## 🚀 Ways to use PANDA

PANDA can be used at three levels. **Installation is optional** for users who
only use the hosted GUI.

1. **Hosted GUI:** open the Hugging Face Space linked above. No local R
   installation is required.
2. **Local GUI:** clone the repository and run `app.R` when local files,
   privacy, or interactive exploration require it.
3. **CLI:** restore the project environment and install the user-level
   `panda` launcher for scripted or server-side analyses.

## 📦 Optional local installation

The following steps are only for local GUI or CLI use. PANDA requires R 4.6.1
and the project dependencies recorded in `renv.lock`.

From a cloned repository:

```bash
git clone https://github.com/kubo-azu/PANDA.git
cd PANDA
Rscript -e 'renv::restore()'
./install_panda.sh
```

The installer places a user-level `panda` launcher in `~/.local/bin`. If that
directory is not on `PATH`, run:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

No administrator privileges are required. The installer does not modify the
system R installation or install packages outside the project environment.

## 🖥️ GUI workflow

From the project directory:

```bash
Rscript -e 'shiny::runApp("app.R")'
```

The GUI supports single-sample and batch workflows for Sanger and amplicon-NGS
data, motif filtering, ASM visualization, heterogeneity statistics, and group
comparison. The GUI is useful for interactive exploration; the CLI is intended
for scripted or server-side analyses.

## 💻 CLI workflow

### 1. Analyze

The usual interface is direct option specification. The minimum NGS command
is:

```bash
panda analyze \
  --mode ngs \
  --input data/fastq \
  --reference data/reference.fasta \
  --output-dir results/my_run
```

The minimum Sanger command is:

```bash
panda analyze \
  --mode sanger \
  --input data/sanger \
  --reference data/reference.fasta \
  --output-dir results/my_sanger_run
```

The main analysis options are:

| Option | Meaning | Default |
|---|---|---|
| `--mode` | `sanger`, `amplicon`, or `ngs` | `amplicon` |
| `--input` | FASTA/FASTQ/AB1 file or input directory | required |
| `--reference` | amplicon reference FASTA | required |
| `--reference-target` | record name when the reference has multiple targets | single reference record |
| `--output-dir` | directory for all result files | `panda_cli_results` |
| `--min-identity` | minimum alignment identity (%) | `90` |
| `--min-conversion` | minimum bisulfite conversion (%) | `95` |
| `--min-count` | minimum NGS multiplicity | `1` |
| `--workers` | parallel alignment workers | `1` |
| `--read-mode` | `merged` or `unmerged` paired-end handling | `merged` |
| `--max-reads` | optional computational read limit | all reads |
| `--max-unique-reads` | optional unique-sequence limit | all unique sequences |
| `--motif-file` | one required DNA motif per line | none |
| `--motifs` | comma-separated required motifs | none |
| `--ab1-trim-start/end` | trim bases from AB1 sequence ends | `20` |
| `--cluster-method` | heterogeneity clustering method | `kmeans` |
| `--k` | number of clusters for k-means | `2` |
| `--seed` | reproducibility seed for clustering | `11` |

For example, a filtered amplicon-NGS run can be written explicitly as:

```bash
panda analyze \
  --mode ngs \
  --input data/fastq \
  --reference data/reference.fasta \
  --reference-target H19_DMR \
  --output-dir results/h19_filtered \
  --min-identity 90 \
  --min-conversion 95 \
  --min-count 2 \
  --workers 8 \
  --motif-file motifs/haplotype_A.txt
```

The most important defaults are:

```text
min_identity       90
min_conversion     95
min_count          1
read_mode          merged
max_reads          all reads
workers            1
```

For long projects or batch reruns, the same options can be stored in a JSON,
YAML, or YML configuration file and passed as one argument. This is optional;
it is a reproducibility convenience, not the normal syntax that users need to
learn first:

```bash
panda analyze --config analysis.yml
```

The demonstration configurations under `examples/` are for reproducing the
paper's synthetic datasets; they are not required to understand or use the
general CLI. Paths are interpreted from the directory in which `panda` is
run, so running it from the repository root is recommended.

For paired-end reads with an unsequenced central region, use
`read_mode: unmerged` and names containing `_R1/_R2` or `_1/_2`. For multiple
reference targets, `reference_target` must be specified.

### 2. Plot

The plotting command accepts any directory containing `PANDA_analysis.json`.
With no figure selection, all figures are created.

```bash
panda plot RESULTS_DIR
```

To select figures, list one or more names after `--figures`:

```bash
panda plot RESULTS_DIR --figures asm samples
```

Available figure names are:

```text
distribution   methylation distribution
heatmap        binned abundance heatmap
lollipop       top-N abundance-ranked methylation patterns
asm            per-sample ASM cluster profile
samples        per-sample overall methylation summary
group          configured group means and SDs
```

`--top-n` is a display option for lollipop plots. It does not change the
analysis or heterogeneity metrics.

### 3. Compare groups

Group comparison is performed after analysis. The recommended input is one
text file per group, with one analyzed sample ID per line. Blank lines and
lines beginning with `#` are ignored.

```bash
panda compare RESULTS_DIR \
  --group-a-file groups/WT.txt \
  --group-b-file groups/KO.txt \
  --name-a WT --name-b KO
```

The command writes CpG-level differences, statistics, and comparison PDFs.
Biological replicates should be independent samples. Comparisons with fewer
than two biological replicates per group are exploratory and should not be
interpreted as replicated statistical inference.

## ⚠️ Input and interpretation requirements

- Use an amplicon reference sequence, not a whole genome.
- Provide the same target and comparable preprocessing across groups.
- Perform adapter/primer trimming and experiment-specific read QC before
  analysis when appropriate.
- Treat technical replicates as technical replicates, not independent
  biological samples.
- Interpret PDR, epipolymorphism, and qFDRP in the context of coverage,
  conversion, identity, and the experimental design.
- Motif filtering is literal sequence filtering; all requested motifs must be
  present in a retained read.

### Sanger input

With `mode: sanger`, PANDA accepts FASTA, multi-FASTA, and AB1-derived input.
Each sequence record in a multi-FASTA file is treated as one clone. When a
directory is supplied, each supported input file is processed as a sample;
multiple FASTA files can therefore be analyzed as a batch. Sanger sequences
are not dereplicated and are not weighted by an NGS-style read count.

For AB1 files, chromatogram base calls are converted to sequences and the
optional `ab1_trim_start` and `ab1_trim_end` settings remove bases at the two
ends before alignment. Primer trimming before analysis remains recommended.

### Amplicon-NGS input

With `mode: amplicon` or `mode: ngs`, PANDA accepts FASTA and FASTQ files,
including `.gz`-compressed FASTQ. FASTQ quality trimming and adapter/primer
removal should normally be performed before PANDA. NGS reads are dereplicated;
the number of original reads represented by each unique sequence is retained
as `Count` and used for abundance-weighted summaries and plots. The default is
to use all reads. `max_reads` and `max_unique_reads` are optional computational
limits, not replacements for experimental QC.

For pre-merged or single-end reads, use `read_mode: merged`. For unmerged
paired-end data, use `read_mode: unmerged`; files must contain `_R1`/`_R2` or
`_1`/`_2` in their names. R2 is reverse-complemented and the pair is reported
as one molecule after the two strands are aligned. Unmerged mode is intended
for long amplicons with an unsequenced central region; overlapping pairs should
normally be merged before analysis to avoid double-counting.

In Sanger mode, `min_count` does not perform dereplication. In NGS modes it is
the minimum multiplicity required for a unique sequence to be retained.

## 🧪 Demonstration data

`PANDA_Demo/` contains synthetic Sanger and amplicon-NGS datasets for ASM,
motif-based haplotype filtering, and a WT/KO group comparison. The archive
`PANDA_Demo.zip` is retained for the GUI download button. Configuration files
under `examples/` are optional reproducibility examples, not analysis results.
See [examples/README.md](examples/README.md) for the recommended examples.

## 🔬 Development and reproducibility

The project uses `renv`; keep `renv.lock`, `renv/activate.R`, `.Rprofile`, and
the source files under version control. Do not commit local `results/`, logs,
RStudio state, or the renv package library.

The analysis implementation is under active development. Before publication,
validate the current release on independent real datasets and report the
alignment, filtering, replicate, and statistical assumptions appropriate to
the experiment.

## 📄 License

See [PANDAcore/LICENSE](PANDAcore/LICENSE).
