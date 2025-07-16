"""Utilities - A collection of Python and shell utilities for app development"""

__version__ = '1.0.4'

# Import all utility modules for easy access (optional imports)
try:
    from . import (
        arg, chat, data, dictionary, error, function, lists, 
        log, multiline, obj, package, parameterize, path, process, 
        python, relation, test
    )
except ImportError as e:
    # Some modules might have optional dependencies
    print(f"Warning: Some modules could not be imported: {e}")

# Common imports (optional)
try:
    from .file import Structure, dups
    from .relation import Relation
    from .lists import nestedly
except ImportError as e:
    # These might fail due to missing dependencies like libmagic
    print(f"Warning: Some utilities could not be imported: {e}")

__all__ = [
    'arg', 'chat', 'data', 'dictionary', 'error', 'function', 'lists',
    'log', 'multiline', 'obj', 'package', 'parameterize', 'path', 'process',
    'python', 'relation', 'test'
]
