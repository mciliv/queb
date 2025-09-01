#!/bin/bash
set -e

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../.. && pwd)"
cd "$ROOT_DIR"

# Actual deploy entrypoint logic goes here.
# Use the main deployment script.

if [ -x scripts/deploy ]; then
  exec scripts/deploy "$@"
fi

echo "No deploy implementation found. Main deployment script should be at scripts/deploy" >&2
exit 1