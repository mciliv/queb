#!/bin/bash

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RC_FILE=~/.zshrc
[[ -f ~/.bashrc ]] && RC_FILE=~/.bashrc

# Remove existing aliases and add new ones
sed -i.bak '/# mol aliases/,/# end mol aliases/d' "$RC_FILE"

cat >> "$RC_FILE" << EOF

# mol aliases
alias dev="$PROJECT_DIR/dev"
alias test="$PROJECT_DIR/test-script"
alias deploy="$PROJECT_DIR/ship"
alias format="$PROJECT_DIR/format"
alias server="$PROJECT_DIR/server"
alias debug="$PROJECT_DIR/debug"
alias unit="$PROJECT_DIR/test-run unit"
alias integration="$PROJECT_DIR/test-run integration"
alias system="$PROJECT_DIR/test-run system"
alias watch="$PROJECT_DIR/test-run watch"
alias pytest="$PROJECT_DIR/test-run pytest"
alias ship="$PROJECT_DIR/ship"
alias ip="$PROJECT_DIR/infrastructure/deployment/ip.sh"
alias mobile="$PROJECT_DIR/infrastructure/deployment/mobile.sh"
alias cert="$PROJECT_DIR/infrastructure/scripts/trust-cert.sh"
alias commit="git commit -am 'Brief description' && git push origin main"
# end mol aliases
EOF

echo "Aliases added to $RC_FILE"
echo "Run 'source $RC_FILE' to reload" 