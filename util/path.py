from pathlib import Path

def prepend(dir_p, ps) -> set[Path]:
    return {Path(dir_p) / p for p in ps}

def top_level(p: Path) -> bool:
    return p == Path(p.root) if p.is_absolute() else not p.parts
