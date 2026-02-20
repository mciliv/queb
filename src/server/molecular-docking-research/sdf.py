#!/usr/bin/env python3
"""Delegates to src/server/sdf.py — see that file for the full implementation."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from sdf import main
if __name__ == "__main__":
    main()
