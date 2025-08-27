#!/bin/bash
set -e

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../.. && pwd)"
cd "$ROOT_DIR"

# Actual deploy entrypoint logic goes here.
# For now, keep it simple and invoke a project-specific helper if present.

if [ -x infrastructure/scripts/project-specific/template-generic.sh ]; then
  exec infrastructure/scripts/project-specific/template-generic.sh "$@"
fi

echo "No deploy implementation found. Add one at infrastructure/scripts/project-specific/template-generic.sh" >&2
exit 1