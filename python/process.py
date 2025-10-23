import time
import os
import logging
import subprocess
from pathlib import Path

from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer


def limit(process: subprocess.Popen, sec_limit: float, output_path: Path):
    stopped = False
    sec = 0
    with output_path.open("a") as output_file:
        while sec < sec_limit:
            time.sleep(1)
            poll = process.poll()
            #python>3.8: if (poll := process.poll()) is not None:
            if poll is not None:
                output_file.write(
                        f"Ran for {sec} sec; process poll value is {poll}"
                        )
                stopped = True
                logging.info(f"Finished: {Path(output_file.name).stem}")
                break
            sec += 1
        if not stopped:
            output_file.write(f"Terminated with {sec_limit} sec limit")


class CustomEventHandler(FileSystemEventHandler):
    def __init__(self, callback):
        super().__init__()
        self.callback = callback

    def on_modified(self, event):
        self.callback(Path(event.src_path))


def call_upon_file_addition(directory_to_watch, callback):
    observer = Observer()
    observer.schedule(CustomEventHandler(callback), str(directory_to_watch), recursive=True)
    observer.start()
    return observer


def recurse_on_subpaths(f: callable[[Path], str], dir_path: Path) -> Path:
    for subpath in dir_path.iterdir(): f(subpath)
    return dirpath

