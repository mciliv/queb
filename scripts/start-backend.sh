#!/bin/bash
set -e

# Start backend with nodemon (shell-based)
# Honors PORT and NODE_ENV env vars if set; provides sensible defaults
PORT="${PORT:-3000}"
NODE_ENV="${NODE_ENV:-development}"

export PORT NODE_ENV

exec npx nodemon \
  --quiet \
  --watch backend \
  --ext js,json \
  --ignore frontend/** \
  --ignore test/** \
  --exec "node --inspect=9229 backend/api/server.js"



