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

The released demonstration datasets are ready to analyze. The optional
`PANDA_demo_data_generation.R` development script requires the mouse BSgenome
package only when regenerating those datasets from scratch; this large
genome-generation dependency is intentionally excluded from the standard
PANDA installation.

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
3. **CLI:** run the installer to restore the project environment and create the
   user-level `panda` command for scripted or server-side analyses.

## 📦 Optional local installation

The following steps are only for local GUI or CLI use. PANDA requires R 4.6.1
and the project dependencies recorded in `renv.lock`. The CLI is supported on
macOS and Linux; Windows users should run it in WSL2.

From a cloned repository, the same installer command is used on macOS, Linux,
and WSL2:

```bash
git clone https://github.com/kubo-azu/PANDA.git
cd PANDA
./install_panda.sh
panda --help
```

The installer checks the selected R version and CPU architecture, restores the
project library, verifies that PANDAcore can be loaded, and creates the
user-level `panda` command. If multiple R installations are available, select
one explicitly:

```bash
./install_panda.sh --rscript /absolute/path/to/Rscript
```

The validated R interpreter is recorded for PANDA only; the installer does not
change the system-wide R selection used by rig, conda, RStudio, or other
projects.

### macOS prerequisites

On macOS, use an R installation whose architecture matches the computer and
the compiled libraries used by R. In particular, Apple Silicon Macs should use
the native arm64 build of R rather than an x86_64 build running through
Rosetta. Check the environment before restoring the project:

```bash
uname -m
Rscript --vanilla -e '
cat(R.version.string, "\n")
cat("R architecture:", R.version$arch, "\n")
cat("R home:", R.home(), "\n")
'
```

For an Apple Silicon installation, `uname -m` should report `arm64`, and R
should report an ARM architecture such as `aarch64` rather than `x86_64`.
Install Apple's Command Line Tools before building any R packages from source:

```bash
xcode-select --install
```

PANDA's R package versions are managed by `renv`, but compilers and any
external libraries needed to build source packages are supplied by the
operating system or an environment manager. The standard PANDA installation
excludes optional demo-generation dependencies, including packages that would
otherwise bring in RCurl and XML. If another package still has to be compiled
from source, an error mentioning a missing library indicates a system
dependency rather than an error in PANDA's methylation analysis.

The native CRAN arm64 R distribution and its matching toolchain are the
recommended setup for Apple Silicon. Do not mix arm64 and x86_64 installations
of R, Homebrew, or compiled R packages. If an older project library was created
with a different R or CPU architecture, restore PANDA in a new project library
rather than reusing the incompatible compiled packages.

The installer places a user-level `panda` launcher in `~/.local/bin`. If that
directory is not on `PATH`, run:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

The installer itself does not require administrator privileges, modify the
system R installation, or install R packages outside the project environment.
Installation of a missing compiler or system library may still require a
system administrator.

### Windows through WSL2

The native Windows CLI is not currently an officially tested target. Windows
users should install WSL2, open its Linux shell, and then follow the same
installation commands shown above. Keep the repository and analysis data in
the WSL Linux filesystem (for example, `~/PANDA` and `~/data`) rather than
under `/mnt/c` when possible, because cross-filesystem access can substantially
reduce I/O performance.

### Optional: use PANDA inside a conda environment

PANDA can also be installed inside a conda environment. This is useful on
servers or when the R runtime itself must be isolated. The recommended
division of responsibility is:

- **conda** manages the R runtime and system-level libraries;
- **renv** manages PANDA's R and Bioconductor packages inside the repository.

For example:

```bash
conda create -n panda-r -c conda-forge r-base=4.6.1
conda activate panda-r

cd PANDA
./install_panda.sh
panda --help
```

Confirm that the intended interpreter is active before restoring the project:

```bash
which Rscript
Rscript -e 'cat(R.version.string, "\n"); cat(R.home(), "\n")'
```

If several R installations are available, select the interpreter while
installing PANDA:

```bash
./install_panda.sh --rscript /absolute/path/to/the/intended/Rscript
```

The `panda` launcher activates this repository's `renv` project even when the
command is run from another directory. User-supplied relative input and output
paths are still resolved from the directory in which the command is invoked.

Do not install PANDA's Bioconductor dependencies separately through both
conda and R; let `renv::restore()` install the versions recorded in
`renv.lock`. The conda environment should use the R version specified by the
lockfile. On macOS, the conda architecture (arm64 or x86_64) should also
match the R and compiled-package ecosystem available on the machine. Hosted
GUI users do not need conda, R, or renv at all.

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
| `--min-shared-cpg` | minimum CpGs shared by a read pair for amplicon qFDRP | `4` |
| `--min-window-coverage` | minimum retained-read/clone coverage for a 4-CpG epipolymorphism window | `2` |
| `--workers` | number of parallel alignment workers (OS-independent; backend selected automatically; maximum 16) | `16` |
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
min_shared_cpg     4
min_window_coverage 2
read_mode          merged
max_reads          all reads
workers            16
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
Summary bars and standard deviations are calculated from sample-level overall
methylation values, and the individual sample values are plotted explicitly.
The overall Wilcoxon test is reported only when each group contains at least
two estimable samples; otherwise its value is `NA`. At each CpG, PANDA uses a
sample-level Welch test when replicate coverage permits it and otherwise
reports the pooled read-level Fisher test as a descriptive fallback. The
`P_Value_Source` column records which test supplied every CpG-level value.
Biological replicates should be independent samples, and pooled read-level
tests must not be interpreted as replicated biological inference.

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

### Heterogeneity metrics and eligibility

PANDA reports three complementary, amplicon-level summaries. Amplicon PDR is
the percentage of eligible reads or clones (at least four observed CpGs) that
contain both methylated and unmethylated calls. Window epipolymorphism is
`1 - sum(p_k^2)` for each consecutive four-CpG window, where `p_k` is the
read-abundance-weighted frequency of pattern `k`; PANDA reports the arithmetic
mean across eligible windows. Amplicon qFDRP is the mean normalized Hamming
distance between retained reads over CpGs observed in both members of a pair.

For dereplicated NGS input, `Count` is the observed number of retained reads
supporting an exact sequence variant. It is used as a read-abundance weight and
must not be interpreted as a UMI-corrected molecule count. The qFDRP calculation
is exactly equivalent to expanding variants by `Count`: pairs of reads from the
same exact variant are included in the denominator and have distance zero.
Sanger records have `Count = 1` and are equally weighted.

Metrics that have no eligible records, windows, or pairs are reported as `NA`,
not zero. Each per-sample heterogeneity CSV records a status, eligible counts,
the weighted number of eligible qFDRP pairs, median shared-CpG count, and the
thresholds used. The defaults are four shared CpGs for qFDRP and coverage of two
for a four-CpG epipolymorphism window; both can be changed explicitly from the
GUI or CLI and are recorded in the analysis bundle.

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
