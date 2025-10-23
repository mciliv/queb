#!/usr/bin/env python3

import sys
from pathlib import Path, PurePath
from IPython import get_ipython, embed


sys.path.append(str(Path.cwd().resolve()))


ipython = get_ipython()

if ipython:
    ipython.magic('load_ext autoreload')
    ipython.magic('autoreload 2')

embed()
