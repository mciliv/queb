# Util - Development Toolkits & Utilities

A comprehensive collection of development toolkits and utilities organized by concept and purpose.

## Structure

```
util/
├── dev-toolkit/         # Generic development & deployment toolkit
│   ├── certs/          # SSL certificate generation
│   ├── scripts/        # Shell scripts for deployment
│   ├── dev.js          # Node.js development server
│   └── README.md       # Development toolkit docs
├── python/             # Python utility modules
│   ├── *.py           # Python utility files
│   ├── tests/         # Python tests
│   └── scripts/       # Python scripts
├── web/               # Web-related utilities
│   ├── capture-screenshot.js
│   ├── css-to-inline.js
│   ├── count_tokens.js
│   └── commit-message-setup.md
├── gcloud-func/       # Google Cloud Functions toolkit
│   ├── deploy.sh      # GCF deployment
│   ├── auth.sh        # Authentication
│   └── README.md      # GCF toolkit docs
├── pyproject.toml     # Python package config
└── README.md          # This file
```

## Development Toolkit

The dev-toolkit provides a complete generic development and deployment system that can be used across multiple projects.

### Quick Start

```bash
# From your project directory (assuming util is at same level)
./scripts/dev serve

# Or directly
node ../util/dev-toolkit/dev.js serve
```

### Features

- **Development Server**: Hot-reload frontend/backend with SSL
- **SSL Certificates**: Automatic trusted certificate generation
- **Deployment**: Google Cloud Functions deployment with domain support
- **Testing**: Unit, integration, and system tests
- **Profiles**: Environment-specific configurations
- **Generic**: Works with any Node.js project

### Project Setup

```javascript
// config/project.js
module.exports = {
  name: 'my-project',
  description: 'My awesome project',
  domain: process.env.DOMAIN_NAME || '',
  functionName: process.env.FUNCTION_NAME || 'my-function',
  pythonVersion: process.env.PYTHON_VERSION || '3.12.9'
};
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

## Web Utilities

Collection of web-related utility scripts:

- **capture-screenshot.js** - Webpage screenshot capture
- **css-to-inline.js** - Convert CSS to inline styles
- **count_tokens.js** - Token counting utilities
- **commit-message-setup.md** - Git commit message guidelines

## Google Cloud Functions Toolkit

Specialized toolkit for Google Cloud Functions:

```bash
# Deploy to GCF
./gcloud-func/deploy.sh

# Setup authentication
./gcloud-func/auth.sh

# Domain management
./gcloud-func/domain.sh
```


## Benefits

- **Well-Organized**: Clear conceptual separation
- **Modular**: Use only what you need
- **Reusable**: Generic tools work across projects
- **Maintained**: Centralized utilities for consistency
- **Documented**: Clear examples and usage guides

## Contributing

When adding new utilities:
1. **Python modules** → `python/` directory
2. **Web utilities** → `web/` directory
3. **Dev tools** → `dev-toolkit/` directory
4. **GCF tools** → `gcloud-func/` directory
5. Update README files with new functionality
6. Add tests when applicable
