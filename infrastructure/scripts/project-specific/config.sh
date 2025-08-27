#!/usr/bin/env bash
set -euo pipefail

# Source shared config from infrastructure/scripts/config.sh (one directory up)
# shellcheck disable=SC1090
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"

# Override PROJECT_ROOT to point to the repository root (parent of infrastructure)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../../.. && pwd)"
export PROJECT_ROOT


