# Scripts Directory

## Overview

This directory contains utility scripts for development, deployment, and infrastructure management. The scripts are organized into two categories:

1. **Project-Specific Scripts** (`project-specific/`) - Scripts tailored to this molecular analysis project
2. **Generic Toolkit** (`.mol-toolkit/`) - Reusable utilities that can be used across multiple projects

## Architecture

### Project-Specific Scripts

Scripts in `project-specific/` are tailored to this specific project and use the generic toolkit utilities:

- `check-domain-status.sh` - Domain health check for queb.space
- `helper.sh` - Python environment setup for this project
- `template-generic.sh` - Example of how to use the toolkit

### Generic Toolkit (`.mol-toolkit/`)

The toolkit provides reusable utilities that can be shared across projects:

- **Core Utilities**: Domain, Google Cloud, Python, and logging functions
- **Templates**: Script templates for common operations
- **Configuration**: Base configuration templates

### Configuration Management

#### `config.sh` - Project-Specific Configuration

This file contains settings specific to the molecular analysis project:

- **Domain settings**: queb.space domain configuration
- **Google Cloud settings**: molecular-analysis function settings
- **Python settings**: Version requirements
- **Paths**: Script and project root directories

#### Benefits of This Architecture:

- **🔧 Reusable**: Toolkit works across different projects
- **📝 Maintainable**: Clear separation of generic and project-specific code
- **🚀 Portable**: Easy to adapt toolkit for new projects
- **🐛 Debuggable**: Consistent logging and error handling

## Best Practices for Context Management

### 1. **Toolkit-First Approach**

- Use generic toolkit functions for common operations
- Only write project-specific logic in project-specific scripts
- Leverage the toolkit's logging and error handling

### 2. **Configuration-First Approach**

- All project-specific values go in `config.sh`
- Scripts should never hardcode domain names, regions, etc.
- Use environment variables with sensible defaults

### 3. **Generic Functions**

- Create reusable functions in the toolkit
- Project scripts should be thin wrappers around toolkit functions
- Functions should accept parameters rather than use global variables

### 4. **Logging and Error Handling**

- Use toolkit's `set_error_handling` for consistent error handling
- All scripts should use toolkit logging functions
- Include timestamps and context in logs

### 5. **Documentation**

- Each script should have a clear purpose statement
- Document required environment variables
- Include usage examples

### 6. **Testing**

- Scripts should be testable in isolation
- Use dry-run modes where possible
- Validate configuration before execution

## Scripts

### Setup and Management Scripts

#### `setup-toolkit.sh` - Initialize Toolkit

**Purpose:** Sets up the mol-toolkit as a git submodule
**Usage:** `./scripts/setup-toolkit.sh`

**What it does:**

- Adds mol-toolkit as a git submodule
- Creates project configuration template
- Sets up project-specific scripts directory
- Makes all toolkit scripts executable
- Commits the setup changes

#### `update-toolkit.sh` - Update Toolkit

**Purpose:** Updates the mol-toolkit to the latest version
**Usage:** `./scripts/update-toolkit.sh`

**What it does:**

- Checks for available updates
- Shows what's changing
- Updates the submodule
- Commits the update

#### `init-project.sh` - Initialize New Project

**Purpose:** Sets up a complete new project with toolkit
**Usage:** `./scripts/init-project.sh`

**What it does:**

- Runs setup-toolkit.sh
- Creates basic project structure
- Generates .gitignore, README.md, package.json
- Creates example project-specific script
- Commits all initialization changes

### Project-Specific Scripts

#### `project-specific/check-domain-status.sh` - Domain Health Check

**Purpose:** Comprehensive domain and DNS status check for queb.space
**Usage:** `./scripts/project-specific/check-domain-status.sh`

**What it checks:**

- DNS zone configuration
- Domain verification status
- Cloud Function deployment
- DNS propagation
- Nameserver configuration

**Uses toolkit functions:** `check_domain_resolution`, `check_cloud_function`, `check_gcloud_dns_zone`

#### `project-specific/helper.sh` - Python Environment Setup

**Purpose:** Sets up Python environment with pyenv and poetry for this project
**Usage:** `source scripts/project-specific/helper.sh`

**What it does:**

- Installs pyenv if needed
- Sets up Python version from config
- Installs poetry
- Configures virtual environment

**Uses toolkit functions:** `setup_python_env`, `install_pyenv`, `install_poetry`

#### `project-specific/template-generic.sh` - Toolkit Usage Example

**Purpose:** Shows how to use the toolkit in project scripts
**Usage:** `./scripts/project-specific/template-generic.sh`

**Demonstrates:**

- Loading toolkit utilities
- Using generic functions
- Consistent logging and error handling

### Deployment Scripts

#### `deploy` - Production Server Deployment

**Command:** `npm run deploy`

**What it does:**

1. ✅ Validates environment variables
2. ✅ Runs pre-deployment tests
3. ✅ Deploys to Google Cloud Functions
4. ✅ Shows deployment URL

**Configuration:**

- Function: `molecular-analysis`
- Runtime: `nodejs20`
- Region: `us-central1`
- Memory: `1GB`
- Timeout: `540s`

#### `ship` - Complete Workflow

**Command:** `npm run ship`

**What it does:**

1. ✅ Stages all changes (`git add .`)
2. ✅ Commits with timestamp (`git commit`)
3. ✅ Runs pre-deployment tests
4. ✅ Deploys to Google Cloud Functions
5. ✅ Pushes to git repository (`git push`)

**Same configuration as `deploy`**

## Benefits

- **📖 Readable** - Clear, documented code instead of long strings
- **🔧 Maintainable** - Easy to modify configuration
- **🐛 Debuggable** - Better error handling and logging
- **⚙️ Configurable** - Centralized configuration object
- **📝 Documented** - Clear comments and step-by-step process

## Environment Variables

**Required:**

- `OPENAI_API_KEY` - OpenAI API key for molecular analysis

**Optional:**

- `GOOGLE_CLOUD_PROJECT` - Google Cloud project ID (for URL display)

## Adapting for New Projects

### Quick Setup (Automated)

1. **Clone your project repository**
2. **Run the initialization script**:
   ```bash
   ./scripts/init-project.sh
   ```
3. **Update the configuration**:
   ```bash
   # Edit scripts/config.sh with your project-specific values
   PROJECT_NAME="your-project-name"
   DOMAIN_NAME="your-domain.com"
   DNS_ZONE_NAME="your-dns-zone"
   FUNCTION_NAME="your-function-name"
   ```
4. **Set up environment variables**:
   ```bash
   export OPENAI_API_KEY="your-api-key"
   export GOOGLE_CLOUD_PROJECT="your-project-id"
   ```
5. **Test the setup**:
   ```bash
   ./scripts/project-specific/example.sh
   ```

### Manual Setup (Alternative)

1. **Set up toolkit**: `./scripts/setup-toolkit.sh`
2. **Create project-specific scripts** using toolkit templates
3. **Update environment variables** in your deployment
4. **Test scripts** to ensure they work with your configuration

### Using Toolkit Templates

1. **Copy template**: `cp .mol-toolkit/templates/domain-check.sh.template scripts/project-specific/your-domain-check.sh`
2. **Update configuration**: Modify the script to load your project's `config.sh`
3. **Customize logic**: Add any project-specific checks or logic
4. **Test**: Run the script to verify it works with your configuration

### Customization Points

- **Domain scripts**: Update DNS zone names and domain mappings
- **Function scripts**: Change function names and regions
- **Python scripts**: Adjust version requirements
- **Logging**: Modify log file paths and formats

### Validation

Run your project-specific domain check script to verify your configuration is working correctly.

### Benefits of This Approach

- **Reusable toolkit**: Generic utilities work across projects
- **Project isolation**: Each project has its own configuration and scripts
- **Easy maintenance**: Update toolkit once, benefits all projects
- **Consistent patterns**: Same logging, error handling, and structure everywhere
