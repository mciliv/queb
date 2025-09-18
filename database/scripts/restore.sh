#!/bin/bash

# restore.sh - PostgreSQL database restore script

set -e

# Load utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../utils.sh"

# Database configuration
DB_NAME="${DB_NAME:-mol_users}"
DB_USER="${DB_USER:-mol_user}"
DB_HOST="${DB_HOST:-localhost}"
DB_PORT="${DB_PORT:-5432}"

# Check if backup file is provided
if [ $# -eq 0 ]; then
    log_error "Usage: $0 <backup_file.sql>"
    log_info "Available backups:"
    ls -la ./backups/backup_*.sql 2>/dev/null || log_info "No backups found in ./backups/"
    exit 1
fi

BACKUP_FILE="$1"

if [ ! -f "$BACKUP_FILE" ]; then
    log_error "Backup file not found: $BACKUP_FILE"
    exit 1
fi

log_warn "⚠️  This will overwrite the current database!"
read -p "Are you sure you want to restore from $BACKUP_FILE? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    log_info "Restore cancelled"
    exit 0
fi

log_info "Starting PostgreSQL restore from: $BACKUP_FILE"

# Terminate active connections to the database
log_info "Terminating active connections..."
psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d postgres -c "
    SELECT pg_terminate_backend(pid)
    FROM pg_stat_activity
    WHERE datname = '$DB_NAME' AND pid <> pg_backend_pid();
" 2>/dev/null || log_warn "Could not terminate connections (may not be necessary)"

# Drop and recreate database
log_info "Recreating database..."
psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d postgres -c "
    DROP DATABASE IF EXISTS $DB_NAME;
    CREATE DATABASE $DB_NAME OWNER $DB_USER;
" 2>/dev/null || log_error "Failed to recreate database"

# Restore from backup
log_info "Restoring from backup..."
PGPASSWORD="$DB_PASSWORD" psql \
    -h "$DB_HOST" \
    -p "$DB_PORT" \
    -U "$DB_USER" \
    -d "$DB_NAME" \
    -f "$BACKUP_FILE"

log_info "✅ Database restore completed successfully"
