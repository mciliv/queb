from __future__ import annotations
import tempfile
import argparse
import magic
import difflib
import shutil
from pathlib import Path
from typing import Any, TypeVar
from hashlib import md5

def temp_file(content=None, dir=None):
    file = tempfile.NamedTemporaryFile(mode='w+t', dir=dir, delete=False)
    if content:
        file.write(content)
        file.flush()
        file.seek(0)
    return file
 

def path(f):
    return Path(f.name)


class File:
    def __init__(self, p: Path=None, content=None, dir=None,
                 create=temp_file, folder_mimetype=None):
        if p:
            self.p = p
        else:
            self.f = create(content, dir)
            self.p = path(self.f)
        self.mimetype = magic.from_file(str(self.p), mime=True) \
            if self.p.is_file() else folder_mimetype

    
    def read(self):
        self.f.seek(0)
        content = self.f.read()
        self.f.seek(0)
        return content

    def recurse(self, f: callable[[Path], Any]):
        if self.p.is_dir():
            for p in self.p.iterdir():
                f(p)


Representation = TypeVar('Representation', str, list)
StructureRepresentation = list[Representation]


class Structure:
    """blueprint & representation expect the top-level
    list then inner lists are directories & strings are files w/ their
    contents
    """
    def __init__(self, blueprint: StructureRepresentation=[], dir=None):
        self.blueprint = blueprint
        self.dir = dir
        self.files = self._files()
    
    def _files(self) -> list[File]:
        return [self._file(content) for content in self.blueprint]

    def _file(self, content: Representation) -> File:
        return Dir(content, self.dir) if isinstance(content, list) \
            else File(content=content, dir=self.dir)

    def clean(self):
        for f in self.files:
            f.structure.clean() if isinstance(f, Dir) else f.f.close()
    
    def representation(self) -> StructureRepresentation:
        return [f.structure.representation() if isinstance(f, Dir) \
                else f.f.read() for f in self.files]
        

class Dir(File):
    def __init__(self, blueprint=[], dir=None,
                 create=tempfile.TemporaryDirectory):
        super().__init__(dir=dir, create=create)
        self.structure = Structure(blueprint, dir)
   

def remove_occurances(string, file_path):
    with open(file_path, "r") as input_file:
        file_contents = input_file.readlines()

    modified_contents = [line.replace(string, "") for line in file_contents]

    with open(file_path, "w") as output_file:
        output_file.writelines(modified_contents)


def remove_spaces_and_line_symbols(name: str):
    for symbol in [" ", "_", "-"]:
        name = name.replace(symbol, "")
    return name


def matches_ignoring_spaces_and_line_symbols(name: str, strings: list[str]):
    return remove_spaces_and_line_symbols(name) in [
        remove_spaces_and_line_symbols(string) for string in strings
    ]


def read_file(file_path):
    with open(file_path) as file:
        return file.readlines()


def eq(bs: list[bytes | None]) -> bool:
    return len(set(b and md5(b).hexdigest() for b in bs)) < 2


def eq_contents(p: Path, candidates: list, candidates_read: callable) \
        -> list:
    with open(p, 'rb') as f:
        return list(filter(
            lambda c: eq([f.read(), candidates_read(c)]), candidates))


def files_are_equivalent(file1_path: Path, file2_path: Path, print_diff=True):
    diff = difflib.ndiff(read_file(file1_path), read_file(file2_path))

    def analyze_line(line):
        return line.startswith("  ")

    return all(analyze_line(line) for line in diff)


def renew(path: Path, mkdirs=True, clean=False, match_as_prefix=False):
    if clean:
        if path.is_file():
            originally_file = True
            path.unlink()
        if path.is_dir():
            shutil.rmtree(path)
        if match_as_prefix:
            for match in path.glob('*'):
                renew(match, clean)
    originally_file = False
    if mkdirs:
        if originally_file:
            path.parent.mkdir(parents=True, exist_ok=True)
        else:
            path.mkdir(parents=True, exist_ok=True)
    return path


def local_data_dir(file_attr: str):
    return renew(Path(file_attr).resolve().parent / "data", mkdirs=True)


def dups(p: Path):
    dups = {}
    for dirpath, dirnames, filenames in p.walk():
        dups[dirpath] = {}
        for f in filenames:
            count = 0
            for f2 in filenames:
                if f == f2:
                    count += 1
            if count > 1:
                dups[dirpath][f] = count
    return dups

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path')
    ds = dups(parser.parse_args().path)
