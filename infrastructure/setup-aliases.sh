#!/bin/bash

# Exit on any error
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Function to setup aliases
setup_mol_aliases() {
    # Validate project root
    if [[ ! -d "$PROJECT_ROOT" ]]; then
        print_error "Project root not found at: $PROJECT_ROOT"
        return 1
    fi

    # Check if scripts exist
    local scripts=("dev" "tests" "ship" "server" "debug" "cleanup")
    for script in "${scripts[@]}"; do
        if [[ ! -f "$PROJECT_ROOT/$script" ]]; then
            print_warning "Script not found: $PROJECT_ROOT/$script"
        fi
    done

    # Determine RC file
    RC_FILE=~/.zshrc
    if [[ -f ~/.bashrc ]] && [[ ! -f ~/.zshrc ]]; then
        RC_FILE=~/.bashrc
        print_warning "Using ~/.bashrc instead of ~/.zshrc"
    fi

    # Create backup
    BACKUP_FILE="${RC_FILE}.bak.$(date +%Y%m%d_%H%M%S)"
    cp "$RC_FILE" "$BACKUP_FILE"
    print_status "Created backup: $BACKUP_FILE"

    # Remove existing mol aliases section
    if grep -q "# mol aliases" "$RC_FILE"; then
        sed -i.bak '/# mol aliases/,/# end mol aliases/d' "$RC_FILE"
        print_status "Removed existing mol aliases"
    fi

    # Add new aliases for direct scripts
    cat >> "$RC_FILE" << EOF

# mol aliases - Direct script aliases
MOL_ROOT="$PROJECT_ROOT"
alias dev="\$MOL_ROOT/dev"
alias test="\$MOL_ROOT/tests"
alias tests="\$MOL_ROOT/tests"
alias ship="\$MOL_ROOT/ship" 
alias server="\$MOL_ROOT/server"
alias debug="\$MOL_ROOT/debug"
alias cleanup="\$MOL_ROOT/cleanup"
alias unit="\$MOL_ROOT/tests unit"
alias integration="\$MOL_ROOT/tests integration"
alias system="\$MOL_ROOT/tests system"
alias watch="\$MOL_ROOT/tests watch"
alias pytest="\$MOL_ROOT/tests pytest"
# end mol aliases
EOF

    print_status "Aliases added to $RC_FILE"
    print_status "Project root: $PROJECT_ROOT"
    print_status "Available aliases: dev, test, tests, ship, server, debug, cleanup, unit, integration, system, watch, pytest"
    print_status "Run 'source $RC_FILE' to reload aliases"
}

# Only run setup if script is executed directly (not sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    setup_mol_aliases
fi 
