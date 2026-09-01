#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Install PANDA and its user-level command launcher.

Usage:
  ./install_panda.sh [--rscript PATH] [--force]

Options:
  --rscript PATH  Use this Rscript executable for PANDA.
  --force         Replace an existing non-PANDA launcher.
  --help, -h      Show this help.

Environment variables:
  R_PANDA_RSCRIPT  Alternative way to select Rscript.
  PANDA_INSTALL_DIR
                    Launcher directory (default: ~/.local/bin, with ~/bin
                    as a fallback when ~/.local/bin is not writable).

The installer validates R and the CPU architecture, restores the project
library from renv.lock, verifies PANDAcore, and records the selected R for the
installed `panda` command. Native Windows users should run the CLI in WSL2.
EOF
}

die() {
  echo "Error: $*" >&2
  exit 1
}

warn() {
  echo "Warning: $*" >&2
}

force=0
rscript_option=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rscript)
      [[ $# -ge 2 ]] || die "--rscript requires a path."
      rscript_option="$2"
      shift 2
      ;;
    --force)
      force=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SOURCE="${SCRIPT_DIR}/bin/panda"
LOCKFILE="${SCRIPT_DIR}/renv.lock"

[[ -x "$SOURCE" ]] || die "executable launcher not found: $SOURCE"
[[ -f "$LOCKFILE" ]] || die "renv.lock was not found: $LOCKFILE"

os_name="$(uname -s 2>/dev/null || echo unknown)"
machine_arch="$(uname -m 2>/dev/null || echo unknown)"
is_wsl=0
wsl_version=""
if [[ "$os_name" == "Linux" ]] && [[ -r /proc/sys/kernel/osrelease ]] && \
   grep -qi microsoft /proc/sys/kernel/osrelease; then
  is_wsl=1
  if grep -qi 'wsl2\|microsoft-standard' /proc/sys/kernel/osrelease; then
    wsl_version="2"
  else
    wsl_version="1"
  fi
fi

case "$os_name" in
  Darwin|Linux) ;;
  MINGW*|MSYS*|CYGWIN*)
    die "native Windows CLI installation is not supported. Use PANDA in WSL2."
    ;;
  *)
    die "unsupported operating system: $os_name (supported: macOS, Linux, WSL2)"
    ;;
esac

if [[ -n "$rscript_option" ]]; then
  rscript_candidate="$rscript_option"
elif [[ -n "${R_PANDA_RSCRIPT:-}" ]]; then
  rscript_candidate="$R_PANDA_RSCRIPT"
else
  rscript_candidate="$(command -v Rscript 2>/dev/null || true)"
fi

[[ -n "$rscript_candidate" ]] || {
  echo "Rscript was not found." >&2
  echo "Install R, or select it explicitly:" >&2
  echo "  ./install_panda.sh --rscript /path/to/Rscript" >&2
  exit 1
}
[[ -x "$rscript_candidate" ]] || die "Rscript is not executable: $rscript_candidate"

lock_r_version="$(
  awk '
    /^[[:space:]]*"R":[[:space:]]*\{/ { in_r = 1; next }
    in_r && /"Version":[[:space:]]*"/ {
      line = $0
      sub(/^.*"Version":[[:space:]]*"/, "", line)
      sub(/".*$/, "", line)
      print line
      exit
    }
  ' "$LOCKFILE"
)"
[[ -n "$lock_r_version" ]] || die "could not read the required R version from renv.lock."

r_info="$(
  "$rscript_candidate" --vanilla -e '
    arch <- R.version$arch
    if (is.null(arch) || !nzchar(arch)) arch <- .Platform$r_arch
    cat(as.character(getRversion()), arch, R.home(), sep = "\t")
  ' 2>/dev/null
)" || die "failed to start Rscript: $rscript_candidate"

IFS=$'\t' read -r r_version r_arch r_home <<< "$r_info"
[[ -n "$r_version" && -n "$r_home" ]] || die "could not determine the selected R configuration."

if [[ "$r_version" != "$lock_r_version" ]]; then
  cat >&2 <<EOF
PANDA requires R ${lock_r_version}, but R ${r_version} was selected.
Selected Rscript: ${rscript_candidate}

Select the intended interpreter explicitly, for example:
  ./install_panda.sh --rscript /path/to/Rscript
EOF
  exit 1
fi

case "${machine_arch}:${r_arch}" in
  arm64:*aarch64*|arm64:*arm64*|aarch64:*aarch64*|aarch64:*arm64*|x86_64:*x86_64*|amd64:*x86_64*) ;;
  *)
    cat >&2 <<EOF
The machine and R architectures do not appear to match.
Machine architecture: ${machine_arch}
R architecture:       ${r_arch}

Use an R build compiled for this machine. Do not mix arm64 and x86_64 R
packages or system libraries.
EOF
    exit 1
    ;;
esac

# Record the interpreter belonging to the validated R installation rather
# than whichever Rscript may later appear first on PATH (for example, after a
# conda environment is activated). This remains specific to PANDA and does not
# change the user's system-wide R selection.
registered_rscript="${r_home}/bin/Rscript"
if [[ ! -x "$registered_rscript" ]]; then
  registered_rscript="$rscript_candidate"
fi

environment_label="$os_name"
if [[ "$is_wsl" -eq 1 ]]; then
  environment_label="WSL${wsl_version}/Linux"
fi

echo "PANDA installation check"
echo "  Environment:    $environment_label"
echo "  Machine:        $machine_arch"
echo "  R version:      $r_version"
echo "  R architecture: $r_arch"
echo "  Rscript:        $registered_rscript"

if [[ "$is_wsl" -eq 1 && "$wsl_version" != "2" ]]; then
  warn "WSL1 was detected; WSL2 is the supported Windows CLI environment."
fi

if [[ "$os_name" == "Darwin" ]] && ! xcode-select -p >/dev/null 2>&1; then
  warn "Apple Command Line Tools were not detected."
  echo "Binary R packages may still install, but source packages require them." >&2
  echo "Install them with: xcode-select --install" >&2
fi

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  echo "  Conda:          ${CONDA_PREFIX}"
fi

if [[ -n "${PANDA_INSTALL_DIR:-}" ]]; then
  INSTALL_DIR="$PANDA_INSTALL_DIR"
else
  preferred_dir="${HOME}/.local/bin"
  if mkdir -p "$preferred_dir" 2>/dev/null && [[ -w "$preferred_dir" ]]; then
    INSTALL_DIR="$preferred_dir"
  else
    INSTALL_DIR="${HOME}/bin"
  fi
fi
TARGET="${INSTALL_DIR}/panda"

if ! mkdir -p "$INSTALL_DIR" 2>/dev/null; then
  cat >&2 <<EOF
Cannot create the installation directory: ${INSTALL_DIR}
Try a user-writable directory, for example:
  PANDA_INSTALL_DIR="\$HOME/bin" ./install_panda.sh
EOF
  exit 1
fi
[[ -w "$INSTALL_DIR" ]] || {
  cat >&2 <<EOF
Installation directory is not writable: ${INSTALL_DIR}
Try a user-writable directory, for example:
  PANDA_INSTALL_DIR="\$HOME/bin" ./install_panda.sh
EOF
  exit 1
}

managed_target=0
if [[ -L "$TARGET" ]]; then
  link_value="$(readlink "$TARGET" 2>/dev/null || true)"
  [[ "$link_value" == "$SOURCE" ]] && managed_target=1
elif [[ -f "$TARGET" ]] && grep -q '^# PANDA launcher generated by install_panda.sh$' "$TARGET" 2>/dev/null; then
  managed_target=1
fi

if [[ -e "$TARGET" || -L "$TARGET" ]]; then
  if [[ "$managed_target" -ne 1 && "$force" -ne 1 ]]; then
    echo "Refusing to overwrite existing file: $TARGET" >&2
    echo "Use --force only if you intend to replace it." >&2
    exit 1
  fi
fi

restore_log="$(mktemp "${TMPDIR:-/tmp}/panda-renv-restore.XXXXXX")"
cleanup() {
  rm -f "$restore_log"
}
trap cleanup EXIT

echo
echo "Restoring the PANDA R environment..."
set +e
RENV_CONFIG_NAMESPACES_CHECK=FALSE PANDA_PROJECT_ROOT="$SCRIPT_DIR" \
  "$registered_rscript" --vanilla -e '
    project <- Sys.getenv("PANDA_PROJECT_ROOT")
    source(file.path(project, "renv", "activate.R"))
    renv::restore(project = project, prompt = FALSE)
    if (!requireNamespace("PANDAcore", quietly = TRUE)) {
      stop("PANDAcore could not be loaded after renv::restore().")
    }
    status <- renv::status(project = project)
    if (!isTRUE(status$synchronized)) {
      stop("The PANDA project library is not synchronized with renv.lock.")
    }
  ' 2>&1 | tee "$restore_log"
restore_status=${PIPESTATUS[0]}
set -e

if [[ "$restore_status" -ne 0 ]]; then
  echo >&2
  echo "PANDA dependency installation did not complete." >&2
  case "$os_name" in
    Darwin)
      cat >&2 <<EOF
If the output mentions a missing system library, install the matching native
library or toolchain described by the package error. Also confirm that R,
Homebrew (if used), and compiled packages share one CPU architecture. See the
macOS prerequisites in README.md.
EOF
      ;;
    Linux)
      cat >&2 <<EOF
The failure may be caused by a missing Linux development library. After renv
has bootstrapped, list the distribution-specific requirements with:
  ${registered_rscript} --vanilla -e 'source("renv/activate.R"); renv::sysreqs(report = TRUE)'
Install system packages through the system administrator, then rerun this
installer.
EOF
      ;;
  esac
  echo "The complete error output was shown above." >&2
  exit "$restore_status"
fi

if [[ -e "$TARGET" || -L "$TARGET" ]]; then
  rm -f "$TARGET"
fi

# Generate a small user-level wrapper. `%q` safely quotes paths containing
# spaces, and the selected interpreter is isolated to PANDA.
{
  echo '#!/usr/bin/env bash'
  echo '# PANDA launcher generated by install_panda.sh'
  echo 'set -euo pipefail'
  printf 'export R_PANDA_RSCRIPT=%q\n' "$registered_rscript"
  printf 'exec %q "$@"\n' "$SOURCE"
} > "$TARGET"
chmod 0755 "$TARGET"

echo
echo "PANDA installation completed."
echo "  Command: $TARGET"
echo "  Rscript: $registered_rscript"

if [[ ":${PATH}:" != *":${INSTALL_DIR}:"* ]]; then
  echo
  echo "Add the command directory to PATH:"
  echo "  export PATH=\"${INSTALL_DIR}:\$PATH\""
  echo "Then run: panda --help"
else
  "$TARGET" --version
  echo "Try: panda --help"
fi
