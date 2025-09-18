#!/bin/bash

# health-check.sh - PostgreSQL database health check script

# Load utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../utils.sh"

# Database configuration
DB_NAME="${DB_NAME:-mol_users}"
DB_USER="${DB_USER:-mol_user}"
DB_HOST="${DB_HOST:-localhost}"
DB_PORT="${DB_PORT:-5432}"

log_info "🔍 Checking PostgreSQL database health..."

# Check connection
if PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -c "SELECT 1;" >/dev/null 2>&1; then
    log_info "✅ Database connection: OK"
else
    log_error "❌ Database connection: FAILED"
    exit 1
fi

# Check database size
DB_SIZE=$(PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -c "
    SELECT pg_size_pretty(pg_database_size('$DB_NAME'));
" 2>/dev/null | tr -d ' ')

if [ -n "$DB_SIZE" ]; then
    log_info "📊 Database size: $DB_SIZE"
else
    log_warn "⚠️  Could not determine database size"
fi

# Check active connections
CONNECTIONS=$(PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -c "
    SELECT COUNT(*) FROM pg_stat_activity WHERE datname = '$DB_NAME';
" 2>/dev/null | tr -d ' ')

if [ -n "$CONNECTIONS" ]; then
    log_info "🔗 Active connections: $CONNECTIONS"
else
    log_warn "⚠️  Could not check connections"
fi

# Check table count
TABLE_COUNT=$(PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -c "
    SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = 'public';
" 2>/dev/null | tr -d ' ')

if [ -n "$TABLE_COUNT" ]; then
    log_info "📋 Tables: $TABLE_COUNT"
else
    log_warn "⚠️  Could not count tables"
fi

# Check recent activity (last 24 hours)
RECENT_ACTIVITY=$(PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -c "
    SELECT COUNT(*) FROM users WHERE last_used > NOW() - INTERVAL '24 hours';
" 2>/dev/null | tr -d ' ')

if [ -n "$RECENT_ACTIVITY" ]; then
    log_info "📈 Recent activity (24h): $RECENT_ACTIVITY users"
else
    log_info "📈 Recent activity: Unable to check"
fi

log_info "🎉 Health check completed!"
