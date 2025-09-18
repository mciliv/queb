#!/bin/bash

# backup.sh - PostgreSQL database backup script

set -e

# Load utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../utils.sh"

# Database configuration
DB_NAME="${DB_NAME:-mol_users}"
DB_USER="${DB_USER:-mol_user}"
DB_HOST="${DB_HOST:-localhost}"
DB_PORT="${DB_PORT:-5432}"

# Backup configuration
BACKUP_DIR="${BACKUP_DIR:-./backups}"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
BACKUP_FILE="$BACKUP_DIR/backup_$TIMESTAMP.sql"

log_info "Starting PostgreSQL backup..."

# Create backup directory if it doesn't exist
mkdir -p "$BACKUP_DIR"

# Perform backup
log_info "Creating backup: $BACKUP_FILE"
PGPASSWORD="$DB_PASSWORD" pg_dump \
    -h "$DB_HOST" \
    -p "$DB_PORT" \
    -U "$DB_USER" \
    -d "$DB_NAME" \
    -f "$BACKUP_FILE"

log_info "✅ Backup completed: $BACKUP_FILE"

# Optional: Clean up old backups (keep last 7 days)
find "$BACKUP_DIR" -name "backup_*.sql" -mtime +7 -delete 2>/dev/null || true
