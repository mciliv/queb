#!/bin/bash

# setup-database.sh - Local PostgreSQL setup for molecular analysis app

set -e

echo "🔧 Setting up PostgreSQL for molecular analysis app..."

# Default values
DB_NAME="mol_users"
DB_USER="mol_user"
DB_PASSWORD="mol_password"
DB_HOST="localhost"
DB_PORT="5432"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

log_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

log_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

log_error() {
    echo -e "${RED}❌ $1${NC}"
}

# Check if PostgreSQL is installed
check_postgresql() {
    log_info "Checking PostgreSQL installation..."
    
    if command -v psql >/dev/null 2>&1; then
        log_success "PostgreSQL client found"
    else
        log_error "PostgreSQL client not found"
        log_info "Please install PostgreSQL:"
        echo "  macOS: brew install postgresql"
        echo "  Ubuntu: sudo apt-get install postgresql postgresql-contrib"
        echo "  Windows: Download from https://www.postgresql.org/download/"
        exit 1
    fi
    
    # Check if PostgreSQL server is running
    if pg_isready -h $DB_HOST -p $DB_PORT >/dev/null 2>&1; then
        log_success "PostgreSQL server is running"
    else
        log_error "PostgreSQL server is not running"
        log_info "Start PostgreSQL server:"
        echo "  macOS: brew services start postgresql"
        echo "  Ubuntu: sudo systemctl start postgresql"
        echo "  Windows: Start from Services or pgAdmin"
        exit 1
    fi
}

# Create database and user
setup_database() {
    log_info "Setting up database and user..."
    
    # Check if running as postgres user or if we have superuser access
    if whoami | grep -q postgres || psql -h $DB_HOST -p $DB_PORT -U postgres -c '\l' >/dev/null 2>&1; then
        POSTGRES_USER="postgres"
    else
        # Try with current user (common on macOS with Homebrew)
        POSTGRES_USER=$(whoami)
        log_info "Attempting to connect as user: $POSTGRES_USER"
    fi
    
    # Create user if it doesn't exist
    log_info "Creating database user: $DB_USER"
    psql -h $DB_HOST -p $DB_PORT -U $POSTGRES_USER -c "
        DO \$\$
        BEGIN
            IF NOT EXISTS (SELECT 1 FROM pg_roles WHERE rolname = '$DB_USER') THEN
                CREATE USER $DB_USER WITH PASSWORD '$DB_PASSWORD';
                ALTER USER $DB_USER CREATEDB;
            END IF;
        END
        \$\$;
    " 2>/dev/null || {
        log_warning "Could not create user with $POSTGRES_USER, trying alternative method..."
        
        # Try creating database directly with current user
        if createdb $DB_NAME 2>/dev/null; then
            log_success "Database created with current user"
            return 0
        else
            log_error "Failed to create database. Please run manually:"
            echo "  sudo -u postgres createuser -s $DB_USER"
            echo "  sudo -u postgres createdb -O $DB_USER $DB_NAME"
            exit 1
        fi
    }
    
    # Create database if it doesn't exist
    log_info "Creating database: $DB_NAME"
    psql -h $DB_HOST -p $DB_PORT -U $POSTGRES_USER -c "
        SELECT 'CREATE DATABASE $DB_NAME OWNER $DB_USER'
        WHERE NOT EXISTS (SELECT FROM pg_database WHERE datname = '$DB_NAME')\gexec
    " 2>/dev/null || {
        log_warning "Could not create database with $POSTGRES_USER, database may already exist"
    }
    
    # Grant privileges
    log_info "Granting privileges..."
    psql -h $DB_HOST -p $DB_PORT -U $POSTGRES_USER -c "
        GRANT ALL PRIVILEGES ON DATABASE $DB_NAME TO $DB_USER;
    " 2>/dev/null || {
        log_warning "Could not grant privileges, may already be set"
    }
}

# Test database connection
test_connection() {
    log_info "Testing database connection..."
    
    if PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -p $DB_PORT -U $DB_USER -d $DB_NAME -c 'SELECT NOW();' >/dev/null 2>&1; then
        log_success "Database connection successful"
    else
        log_error "Database connection failed"
        log_info "Please check your PostgreSQL installation and try again"
        exit 1
    fi
}

# Create .env file for local development
create_env_file() {
    log_info "Creating .env file for local development..."
    
    ENV_FILE="$(dirname "$(dirname "$(pwd)")")"/.env
    
    if [ ! -f "$ENV_FILE" ]; then
        cat > "$ENV_FILE" << EOF
# PostgreSQL Configuration for Local Development
DB_HOST=$DB_HOST
DB_PORT=$DB_PORT
DB_NAME=$DB_NAME
DB_USER=$DB_USER
DB_PASSWORD=$DB_PASSWORD

# Development settings
NODE_ENV=development
OPENAI_API_KEY=your_openai_api_key_here
STRIPE_PUBLISHABLE_KEY=your_stripe_publishable_key_here
STRIPE_SECRET_KEY=your_stripe_secret_key_here
EOF
        log_success ".env file created at $ENV_FILE"
        log_warning "Please update the API keys in .env file"
    else
        log_info ".env file already exists, skipping creation"
    fi
}

# Main execution
main() {
    echo "🧬 Molecular Analysis App - Database Setup"
    echo "=========================================="
    
    check_postgresql
    setup_database
    test_connection
    create_env_file
    
    echo ""
    log_success "Database setup completed successfully!"
    echo ""
    log_info "Connection details:"
    echo "  Host: $DB_HOST"
    echo "  Port: $DB_PORT"
    echo "  Database: $DB_NAME"
    echo "  User: $DB_USER"
    echo ""
    log_info "To connect manually:"
    echo "  PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -p $DB_PORT -U $DB_USER -d $DB_NAME"
    echo ""
    log_info "Your app will automatically create the required tables on startup."
}

# Run main function
main "$@" 