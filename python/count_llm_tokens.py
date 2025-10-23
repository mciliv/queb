#!/usr/bin/env python3
import argparse
import dataclasses
import datetime
import fnmatch
import json
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


def resolve_encoding(encoding_name: str):
    import tiktoken

    known_encodings = {
        "cl100k_base",
        "o200k_base",
        "p50k_base",
        "r50k_base",
    }
    if encoding_name not in known_encodings:
        raise SystemExit(
            f"Unsupported encoding '{encoding_name}'. Choose one of: "
            + ", ".join(sorted(known_encodings))
        )
    return tiktoken.get_encoding(encoding_name)


@dataclasses.dataclass
class RepoStats:
    name: str
    path: Path
    files_counted: int = 0
    bytes_counted: int = 0
    tokens: int = 0


DEFAULT_EXCLUDE_DIR_NAMES: Set[str] = {
    ".git",
    "node_modules",
    "dist",
    "build",
    ".next",
    ".venv",
    "venv",
    "__pycache__",
    ".cache",
    ".idea",
    ".vscode",
    ".pytest_cache",
    "coverage",
    "logs",
    "log",
    "tmp",
    "out",
}

DEFAULT_EXCLUDE_FILE_GLOBS: Sequence[str] = (
    "*.jpg",
    "*.jpeg",
    "*.png",
    "*.gif",
    "*.webp",
    "*.svg",
    "*.ico",
    "*.pdf",
    "*.zip",
    "*.tar",
    "*.gz",
    "*.rar",
    "*.7z",
    "*.mp3",
    "*.ogg",
    "*.wav",
    "*.mp4",
    "*.mov",
    "*.avi",
    "*.mkv",
    "*.exe",
    "*.dll",
    "*.so",
    "*.dylib",
    "*.otf",
    "*.ttf",
    "*.woff",
    "*.woff2",
    "*.eot",
    "*.db",
    "*.sqlite",
    "*.sqlite3",
    "*.pyc",
    "*.pkl",
    "*.bin",
)


def is_hidden_path_component(name: str) -> bool:
    return name.startswith(".")


def file_is_excluded(file_path: Path, exclude_file_globs: Sequence[str]) -> bool:
    file_name = file_path.name
    for pattern in exclude_file_globs:
        if fnmatch.fnmatch(file_name, pattern):
            return True
    return False


def iter_repo_dirs(root: Path) -> List[Path]:
    repo_dirs: List[Path] = []
    for child in root.iterdir():
        if not child.is_dir():
            continue
        # Treat immediate children that are git repos as repos
        if (child / ".git").is_dir():
            repo_dirs.append(child)
    return repo_dirs


def count_tokens_in_dir(
    dir_path: Path,
    encoding_name: str,
    exclude_dir_names: Set[str],
    exclude_file_globs: Sequence[str],
    include_hidden: bool,
    max_bytes: int,
    skip_immediate_subdirs: Optional[Set[Path]] = None,
) -> Tuple[int, int, int]:
    encoding = resolve_encoding(encoding_name)
    total_tokens = 0
    total_files = 0
    total_bytes = 0

    dir_path = dir_path.resolve()
    skip_subdir_names: Set[str] = set()
    if skip_immediate_subdirs:
        skip_subdir_names = {p.name for p in skip_immediate_subdirs}

    for current_root, dirnames, filenames in os.walk(dir_path, topdown=True):
        current_root_path = Path(current_root)

        if current_root_path == dir_path and skip_subdir_names:
            dirnames[:] = [d for d in dirnames if d not in skip_subdir_names]

        # Filter directories
        filtered_dirnames: List[str] = []
        for d in dirnames:
            if d in exclude_dir_names:
                continue
            if not include_hidden and is_hidden_path_component(d):
                # Still allow .git to be excluded by default above
                continue
            filtered_dirnames.append(d)
        dirnames[:] = filtered_dirnames

        for fname in filenames:
            if not include_hidden and is_hidden_path_component(fname):
                continue
            file_path = current_root_path / fname
            if file_is_excluded(file_path, exclude_file_globs):
                continue
            try:
                try:
                    size_bytes = file_path.stat().st_size
                except FileNotFoundError:
                    # Skip transient files
                    continue
                if size_bytes > max_bytes:
                    continue
                # Read as text; replace invalid bytes
                with open(file_path, "rb") as f:
                    raw = f.read()
                try:
                    text = raw.decode("utf-8", errors="replace")
                except Exception:
                    # Fallback to latin-1 to cover odd encodings
                    text = raw.decode("latin-1", errors="replace")

                token_ids = encoding.encode(text)
                n_tokens = len(token_ids)
                total_tokens += n_tokens
                total_files += 1
                total_bytes += size_bytes
            except Exception as e:
                print(f"WARN: failed to process {file_path}: {e}", file=sys.stderr)
                continue

    return total_tokens, total_files, total_bytes


def format_int(n: int) -> str:
    return f"{n:,}"


def print_table(rows: List[Tuple[str, int, int]]):
    # rows: (name, files, tokens)
    name_width = max((len(name) for name, _, _ in rows), default=4)
    files_width = max((len(format_int(files)) for _, files, _ in rows), default=5)
    tokens_width = max((len(format_int(tokens)) for _, _, tokens in rows), default=6)
    header = (
        f"{'Repo':{name_width}}  {'Files':>{files_width}}  {'Tokens':>{tokens_width}}"
    )
    sep = f"{'-' * name_width}  {'-' * files_width}  {'-' * tokens_width}"
    print(header)
    print(sep)
    for name, files, tokens in rows:
        print(
            f"{name:{name_width}}  {format_int(files):>{files_width}}  {format_int(tokens):>{tokens_width}}"
        )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Count LLM tokens in a directory (per immediate child git repo and overall)."
        )
    )
    parser.add_argument(
        "--root",
        default=str(Path("~").expanduser() / "Code"),
        help="Root directory to scan (default: ~/Code)",
    )
    parser.add_argument(
        "--encoding",
        default="cl100k_base",
        help="Tokenizer encoding name (cl100k_base, o200k_base, p50k_base, r50k_base)",
    )
    parser.add_argument(
        "--max-bytes",
        type=int,
        default=5 * 1024 * 1024,
        help="Skip files larger than this many bytes (default: 5 MiB)",
    )
    parser.add_argument(
        "--include-hidden",
        action="store_true",
        help="Include hidden files and directories (default: false)",
    )
    parser.add_argument(
        "--extra-exclude-dir",
        action="append",
        default=[],
        help="Additional directory names to exclude (can be repeated)",
    )
    parser.add_argument(
        "--json-out",
        default=None,
        help="Write detailed results to this JSON file",
    )
    parser.add_argument(
        "--csv-out",
        default=None,
        help="Write summary table to this CSV file",
    )

    args = parser.parse_args(argv)

    root = Path(args.root).expanduser().resolve()
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"Root does not exist or is not a directory: {root}")

    # Resolve encoding early to validate
    _ = resolve_encoding(args.encoding)

    exclude_dir_names = set(DEFAULT_EXCLUDE_DIR_NAMES)
    if args.extra_exclude_dir:
        exclude_dir_names.update(args.extra_exclude_dir)

    repo_dirs = iter_repo_dirs(root)
    repo_names = [p.name for p in repo_dirs]

    results: List[RepoStats] = []

    # Compute stats for each repo
    for repo_dir in sorted(repo_dirs, key=lambda p: p.name.lower()):
        tokens, files, bytes_count = count_tokens_in_dir(
            dir_path=repo_dir,
            encoding_name=args.encoding,
            exclude_dir_names=exclude_dir_names,
            exclude_file_globs=DEFAULT_EXCLUDE_FILE_GLOBS,
            include_hidden=args.include_hidden,
            max_bytes=args.max_bytes,
            skip_immediate_subdirs=None,
        )
        results.append(
            RepoStats(
                name=repo_dir.name,
                path=repo_dir,
                files_counted=files,
                bytes_counted=bytes_count,
                tokens=tokens,
            )
        )

    # Count root-level files and non-repo subdirectories as "root"
    skip_subdirs = set(repo_dirs)
    root_tokens, root_files, root_bytes = count_tokens_in_dir(
        dir_path=root,
        encoding_name=args.encoding,
        exclude_dir_names=exclude_dir_names,
        exclude_file_globs=DEFAULT_EXCLUDE_FILE_GLOBS,
        include_hidden=args.include_hidden,
        max_bytes=args.max_bytes,
        skip_immediate_subdirs=skip_subdirs,
    )
    # Name the root bucket clearly
    results.append(
        RepoStats(
            name=root.name,
            path=root,
            files_counted=root_files,
            bytes_counted=root_bytes,
            tokens=root_tokens,
        )
    )

    # Prepare printable table
    summary_rows: List[Tuple[str, int, int]] = [
        (r.name, r.files_counted, r.tokens) for r in results
    ]
    summary_rows.sort(key=lambda t: t[2], reverse=True)

    print(f"Root: {root}")
    print(f"Encoding: {args.encoding}")
    print()
    print_table(summary_rows)

    grand_total_tokens = sum(r.tokens for r in results)
    grand_total_files = sum(r.files_counted for r in results)
    print()
    print(f"Total files: {format_int(grand_total_files)}")
    print(f"Total tokens: {format_int(grand_total_tokens)}")

    # Optional outputs
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    if args.json_out:
        report = {
            "root": str(root),
            "encoding": args.encoding,
            "generated_at": timestamp,
            "repos": [
                {
                    "name": r.name,
                    "path": str(r.path),
                    "files": r.files_counted,
                    "bytes": r.bytes_counted,
                    "tokens": r.tokens,
                }
                for r in results
            ],
            "totals": {
                "files": grand_total_files,
                "tokens": grand_total_tokens,
            },
        }
        with open(args.json_out, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)

    if args.csv_out:
        import csv

        with open(args.csv_out, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["repo", "files", "tokens"])
            for name, files, tokens in summary_rows:
                writer.writerow([name, files, tokens])
            writer.writerow(["TOTAL", grand_total_files, grand_total_tokens])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


