#!/usr/bin/env sh

set -eu
# shellcheck shell=sh

# Best-effort pipefail for bash/zsh; ignored if unsupported in POSIX sh
{ set -o pipefail 2>/dev/null || true; } >/dev/null 2>&1 || true

ROOT_DEFAULT="$HOME/Code"
ROOT_DIR="${1:-$ROOT_DEFAULT}"

VENV_DIR="$ROOT_DIR/.venv-token-count"
SCRIPT_PATH="$ROOT_DIR/util/python/count_llm_tokens.py"

if [ ! -d "$ROOT_DIR" ]; then
  echo "Root directory does not exist: $ROOT_DIR" >&2
  exit 1
fi

if [ ! -f "$SCRIPT_PATH" ]; then
  echo "Token counter script not found: $SCRIPT_PATH" >&2
  exit 1
fi

PYTHON_BIN="python3"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "python3 is required but not found in PATH" >&2
  exit 1
fi

# Create venv if needed
if [ ! -d "$VENV_DIR" ]; then
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

# Activate venv (POSIX compatible)
. "$VENV_DIR/bin/activate"

# Ensure pip is present
python -m ensurepip --upgrade >/dev/null 2>&1 || true
python -m pip install --upgrade pip >/dev/null 2>&1

# Install tiktoken if missing
if ! python -c 'import tiktoken' >/dev/null 2>&1; then
  python -m pip install 'tiktoken==0.7.0' >/dev/null 2>&1
fi

# Outputs
JSON_OUT="$ROOT_DIR/token_counts.json"
CSV_OUT="$ROOT_DIR/token_counts.csv"

# Pass through any additional flags after the first arg (root). To pass a custom root, use:
#   ./count_llm_tokens.sh /path/to/root -- --encoding o200k_base
EXTRA_ARGS=""
if [ "$#" -ge 2 ]; then
  shift 1
  if [ "${1:-}" = "--" ]; then
    shift 1
    EXTRA_ARGS="$*"
  fi
fi

echo "Counting tokens under: $ROOT_DIR"
python "$SCRIPT_PATH" --root "$ROOT_DIR" --encoding cl100k_base \
  --json-out "$JSON_OUT" --csv-out "$CSV_OUT" ${EXTRA_ARGS}

echo
echo "Wrote JSON: $JSON_OUT"
echo "Wrote CSV:  $CSV_OUT"


