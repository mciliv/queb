#!/bin/bash

# Essential Utilities for Development Scripts
# ==========================================


# Simple logging with colors
log_info()  { echo -e "\033[32m✅ $1\033[0m"; }
log_warn()  { echo -e "\033[33m⚠️  $1\033[0m"; }
log_error() { echo -e "\033[31m❌ $1\033[0m"; }

# Essential utility functions
check_port() {
    local port="$1"
    lsof -i ":$port" > /dev/null 2>&1
}

kill_port() {
    local port="$1"
    local pid=$(lsof -ti ":$port" 2>/dev/null)
    if [ -n "$pid" ]; then
        kill -9 "$pid" 2>/dev/null
        return 0
    fi
    return 1
}

check_url() {
    local url="$1"
    curl -s --max-time 5 "$url" > /dev/null 2>&1
} 