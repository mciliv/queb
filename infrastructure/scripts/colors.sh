#!/bin/bash

# Shared color utility for all scripts
# Usage: source infrastructure/scripts/colors.sh

# Color definitions
export RED='\033[0;31m'
export GREEN='\033[0;32m'
export YELLOW='\033[1;33m'
export BLUE='\033[0;34m'
export MAGENTA='\033[0;35m'
export CYAN='\033[0;36m'
export WHITE='\033[1;37m'
export GRAY='\033[0;90m'
export NC='\033[0m'  # No Color

# Logging functions with consistent emojis and colors
log_info() { echo -e "${BLUE}ℹ️  $1${NC}"; }
log_success() { echo -e "${GREEN}✅ $1${NC}"; }
log_warning() { echo -e "${YELLOW}⚠️  $1${NC}"; }
log_error() { echo -e "${RED}❌ $1${NC}"; }
log_debug() { echo -e "${GRAY}🐛 $1${NC}"; }
log_ship() { echo -e "${CYAN}🚢 $1${NC}"; }
log_test() { echo -e "${MAGENTA}🧪 $1${NC}"; }
log_clean() { echo -e "${YELLOW}🧹 $1${NC}"; }
log_rocket() { echo -e "${GREEN}🚀 $1${NC}"; }

# Header function for script titles
log_header() {
    echo -e "${WHITE}======================================${NC}"
    echo -e "${WHITE}$1${NC}"
    echo -e "${WHITE}======================================${NC}"
}

# Simple colored echo without emojis (for when you just want color)
echo_red() { echo -e "${RED}$1${NC}"; }
echo_green() { echo -e "${GREEN}$1${NC}"; }
echo_yellow() { echo -e "${YELLOW}$1${NC}"; }
echo_blue() { echo -e "${BLUE}$1${NC}"; }
echo_magenta() { echo -e "${MAGENTA}$1${NC}"; }
echo_cyan() { echo -e "${CYAN}$1${NC}"; } 