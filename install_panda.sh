#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Install the PANDA command for the current user.

Usage:
  ./install_panda.sh [--force]

The command is installed as ~/.local/bin/panda.
No R packages are installed by this script; use renv::restore() separately.
EOF
}

force=0
case "${1:-}" in
  "") ;;
  --force) force=1 ;;
  --help|-h) usage; exit 0 ;;
  *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
esac

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript was not found on PATH." >&2
  echo "Install R and make sure the intended Rscript is available." >&2
  exit 1
fi

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
SOURCE="${SCRIPT_DIR}/bin/panda"
if [[ ! -x "$SOURCE" ]]; then
  echo "Error: executable launcher not found: $SOURCE" >&2
  exit 1
fi

INSTALL_DIR="${PANDA_INSTALL_DIR:-${HOME}/.local/bin}"
TARGET="${INSTALL_DIR}/panda"
mkdir -p "$INSTALL_DIR"

if [[ -e "$TARGET" || -L "$TARGET" ]]; then
  if [[ "$force" -ne 1 ]]; then
    echo "Refusing to overwrite existing file: $TARGET" >&2
    echo "Use --force only if you intend to replace it." >&2
    exit 1
  fi
  rm -f "$TARGET"
fi

ln -s "$SOURCE" "$TARGET"

echo "Installed PANDA command: $TARGET"
if [[ ":${PATH}:" != *":${INSTALL_DIR}:"* ]]; then
  echo "Add this directory to PATH, for example:"
  echo "  export PATH=\"${INSTALL_DIR}:\$PATH\""
else
  echo "Try: panda --help"
fi
echo "Rscript: $(command -v Rscript)"
