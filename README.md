# Util - Python & Web Development Utilities

A comprehensive collection of Python utilities and a complete webkit for app/web development.

## Structure

```
util/
├── webkit/              # Complete web/app development toolkit
│   ├── deploy/          # Deployment utilities (shell scripts)
│   ├── frontend/        # Frontend utilities (CSS, JS, HTML)
│   ├── backend/         # Backend utilities (Node.js, Python)
│   ├── web/             # Web-specific utilities (domains, SSL, CDN)
│   ├── app/             # App development utilities (mobile, desktop, PWA)
│   └── README.md        # Webkit documentation
├── util/                # Python utility modules
│   ├── *.py            # Python utility files
│   ├── tests/          # Python tests
│   └── scripts/        # Python scripts
├── pyproject.toml       # Python package config
└── README.md           # This file
```

## Python Utilities

### Installation

```bash
# From git
pip install git+https://github.com/mciliv/util.git

# Or add to pyproject.toml
dependencies = [
    "util @ git+https://github.com/mciliv/util.git"
]
```

### Usage

```python
# Import specific utilities
from util.file import Structure, dups
from util.relation import Relation
from util.lists import nestedly

# Or import the whole package
import util
```

### Available Modules

- **file.py** - File operations, Structure class, duplicate detection
- **relation.py** - Relationship management, Relation class
- **lists.py** - List utilities, nested operations
- **chat.py** - Chat/communication utilities
- **dictionary.py** - Dictionary operations
- **process.py** - Process management
- **package.py** - Package utilities
- **log.py** - Logging utilities
- **error.py** - Error handling
- **path.py** - Path utilities
- **python.py** - Python-specific utilities
- **test.py** - Testing utilities
- **arg.py** - Argument parsing
- **data.py** - Data utilities
- **function.py** - Function utilities
- **multiline.py** - Multiline text processing
- **obj.py** - Object utilities
- **parameterize.py** - Parameter utilities

## WebKit - Complete App & Web Development Toolkit

The webkit provides everything needed for modern web and app development:

### Quick Start

```bash
# Copy webkit to your project
cp -r util/webkit/ your-project/

# Use deployment utilities
source webkit/deploy/utils.sh
source webkit/deploy/config.sh

# Deploy your app
./webkit/deploy/deploy.sh
```

### What's Included

#### Deployment
- Domain & DNS management
- Google Cloud deployment
- SSL certificate setup
- CDN configuration

#### Frontend
- CSS utilities and responsive design
- JavaScript utilities (DOM, API, Storage, Validation)
- HTML templates with SEO optimization
- Asset management

#### Backend
- Node.js utilities (Express, Database, Validation, Security)
- Python backend utilities
- API development tools

#### Web
- Domain registration and management
- SSL certificate utilities
- CDN setup and optimization
- Web monitoring and analytics

#### App Development
- Mobile app utilities
- Desktop app utilities
- Progressive Web App (PWA) setup

### Usage Examples

#### Deploy a Web App
```bash
#!/bin/bash
source webkit/deploy/utils.sh
source webkit/deploy/config.sh

setup_node_project
install_node_deps
build_node_app
deploy_function "$FUNCTION_NAME" "$REGION" "$SOURCE_DIR"
setup_domain "$DOMAIN_NAME"
setup_ssl "$DOMAIN_NAME"
```

#### Create a PWA
```bash
#!/bin/bash
source webkit/app/pwa/utils.sh

setup_pwa "$PROJECT_NAME"
generate_pwa_manifest
setup_service_worker
build_node_app
deploy_function "$FUNCTION_NAME" "$REGION" "$SOURCE_DIR"
```

## Benefits

- **Comprehensive**: Covers all aspects of development
- **Modular**: Use only what you need
- **Portable**: Easy to copy between projects
- **Well-documented**: Clear examples and usage
- **Maintained**: Centralized utilities for consistency

## Contributing

When adding new utilities:
1. Add Python modules to `util/` directory
2. Add web/app utilities to appropriate `webkit/` subdirectory
3. Include documentation and examples
4. Update README files with new functionality
5. Add tests if applicable
