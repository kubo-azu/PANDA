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

if [[ -n "${PANDA_INSTALL_DIR:-}" ]]; then
  INSTALL_DIR="${PANDA_INSTALL_DIR}"
else
  # Prefer the conventional per-user location.  Some managed macOS
  # installations contain a root-owned ~/.local/bin; in that case fall back
  # automatically to ~/bin instead of failing with a cryptic ln error.
  PREFERRED_DIR="${HOME}/.local/bin"
  if mkdir -p "$PREFERRED_DIR" 2>/dev/null && [[ -w "$PREFERRED_DIR" ]]; then
    INSTALL_DIR="$PREFERRED_DIR"
  else
    INSTALL_DIR="${HOME}/bin"
  fi
fi
TARGET="${INSTALL_DIR}/panda"
if ! mkdir -p "$INSTALL_DIR" 2>/dev/null; then
  echo "Error: cannot create installation directory: $INSTALL_DIR" >&2
  echo "Try a user-writable directory, for example:" >&2
  echo "  PANDA_INSTALL_DIR=\"\$HOME/bin\" ./install_panda.sh" >&2
  exit 1
fi
if [[ ! -w "$INSTALL_DIR" ]]; then
  echo "Error: installation directory is not writable: $INSTALL_DIR" >&2
  echo "Try a user-writable directory, for example:" >&2
  echo "  PANDA_INSTALL_DIR=\"\$HOME/bin\" ./install_panda.sh" >&2
  exit 1
fi

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
  echo "Then run: panda --help"
else
  echo "Try: panda --help"
fi
echo "Rscript: $(command -v Rscript)"
