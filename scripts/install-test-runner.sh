#!/bin/bash

# Mol Test Runner Installation Script
# Supports: PM2, systemd, Docker

set -e

INSTALL_DIR="/opt/mol"
SERVICE_NAME="mol-test-runner"
LOG_DIR="/var/log/mol"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

check_root() {
    if [[ $EUID -ne 0 ]]; then
        error "This script must be run as root"
        exit 1
    fi
}

detect_os() {
    if [[ -f /etc/os-release ]]; then
        . /etc/os-release
        OS=$ID
        VER=$VERSION_ID
    else
        error "Cannot detect OS"
        exit 1
    fi
}

install_dependencies() {
    log "Installing system dependencies..."
    
    case $OS in
        ubuntu|debian)
            apt-get update
            apt-get install -y nodejs npm python3 python3-pip git curl
            pip3 install rdkit
            ;;
        centos|rhel|fedora)
            yum install -y nodejs npm python3 python3-pip git curl
            pip3 install rdkit
            ;;
        *)
            error "Unsupported OS: $OS"
            exit 1
            ;;
    esac
    
    # Install PM2 globally
    npm install -g pm2
    
    # Install node-cron for scheduler
    npm install -g node-cron
}

setup_directories() {
    log "Setting up directories..."
    
    # Create installation directory
    mkdir -p $INSTALL_DIR
    mkdir -p $LOG_DIR
    mkdir -p $INSTALL_DIR/logs
    mkdir -p $INSTALL_DIR/sdf_files
    
    # Create node user if doesn't exist
    if ! id -u node >/dev/null 2>&1; then
        useradd -r -s /bin/false node
    fi
}

install_application() {
    log "Installing application files..."
    
    # Copy application files
    cp -r ./* $INSTALL_DIR/
    
    # Install dependencies
    cd $INSTALL_DIR
    npm ci --only=production
    
    # Set permissions
    chown -R node:node $INSTALL_DIR
    chown -R node:node $LOG_DIR
}

setup_pm2() {
    log "Setting up PM2..."
    
    # Start PM2 as node user
    sudo -u node PM2_HOME=/home/node/.pm2 pm2 start ecosystem.config.js
    
    # Save PM2 process list
    sudo -u node PM2_HOME=/home/node/.pm2 pm2 save
    
    # Setup PM2 startup script
    pm2 startup systemd -u node --hp /home/node
}

setup_systemd() {
    log "Setting up systemd service..."
    
    # Copy service file
    cp scripts/mol-test-runner.service /etc/systemd/system/
    
    # Update paths in service file
    sed -i "s|/opt/mol|$INSTALL_DIR|g" /etc/systemd/system/mol-test-runner.service
    
    # Reload systemd and enable service
    systemctl daemon-reload
    systemctl enable mol-test-runner.service
    systemctl start mol-test-runner.service
}

setup_docker() {
    log "Setting up Docker container..."
    
    # Check if Docker is installed
    if ! command -v docker &> /dev/null; then
        error "Docker is not installed. Please install Docker first."
        exit 1
    fi
    
    # Build and start containers
    docker-compose -f docker-compose.test.yml up -d --build
}

setup_cron_backup() {
    log "Setting up cron backup..."
    
    # Add cron job for backup test runner
    cat > /etc/cron.d/mol-test-runner << EOF
# Mol Test Runner Cron Backup
# Run smoke tests every 30 minutes as fallback
*/30 * * * * node $INSTALL_DIR/node_modules/.bin/jest --selectProjects smoke >> $LOG_DIR/cron-smoke.log 2>&1

# Daily full test at 3 AM
0 3 * * * node $INSTALL_DIR/scripts/background-test-runner.js run >> $LOG_DIR/cron-daily.log 2>&1

# Health check every 5 minutes
*/5 * * * * node $INSTALL_DIR/scripts/health-check.js || (echo "Health check failed" | mail -s "Mol Test Runner Alert" admin@example.com)
EOF
}

setup_logrotate() {
    log "Setting up log rotation..."
    
    cat > /etc/logrotate.d/mol-test-runner << EOF
$LOG_DIR/*.log {
    daily
    rotate 7
    compress
    delaycompress
    missingok
    notifempty
    create 0644 node node
    sharedscripts
    postrotate
        pm2 reloadLogs
    endscript
}
EOF
}

print_status() {
    log "Installation complete!"
    echo ""
    echo "Test Runner Status:"
    echo "=================="
    
    if command -v pm2 &> /dev/null; then
        sudo -u node PM2_HOME=/home/node/.pm2 pm2 status
    fi
    
    if systemctl is-active --quiet mol-test-runner; then
        echo "Systemd service: Active ✅"
    else
        echo "Systemd service: Inactive ❌"
    fi
    
    if docker ps | grep -q mol-test-runner; then
        echo "Docker container: Running ✅"
    else
        echo "Docker container: Not running ❌"
    fi
    
    echo ""
    echo "Useful commands:"
    echo "- pm2 status                    # Check PM2 processes"
    echo "- systemctl status mol-test-runner  # Check systemd service"
    echo "- docker-compose -f docker-compose.test.yml logs  # Check Docker logs"
    echo "- tail -f $LOG_DIR/*.log        # View logs"
}

# Main installation flow
main() {
    log "Starting Mol Test Runner installation..."
    
    check_root
    detect_os
    
    echo "Choose installation method:"
    echo "1) PM2 (recommended for development)"
    echo "2) systemd (recommended for production)"
    echo "3) Docker (recommended for isolation)"
    echo "4) All methods"
    read -p "Enter choice [1-4]: " choice
    
    install_dependencies
    setup_directories
    install_application
    
    case $choice in
        1)
            setup_pm2
            ;;
        2)
            setup_systemd
            ;;
        3)
            setup_docker
            ;;
        4)
            setup_pm2
            setup_systemd
            setup_docker
            ;;
        *)
            error "Invalid choice"
            exit 1
            ;;
    esac
    
    setup_cron_backup
    setup_logrotate
    print_status
}

# Run main function
main "$@"