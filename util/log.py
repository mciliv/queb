from pathlib import Path
import logging


def setup(dir_path: Path, file_name='logs.log'):
    logging.basicConfig(level=logging.INFO)
    file_handler = logging.FileHandler(str(dir_path / file_name))
    logging.getLogger().addHandler(file_handler)

