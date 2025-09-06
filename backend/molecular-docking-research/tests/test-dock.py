from pathlib import Path

from util.file import dups

def test_dups():
    p = Path(__file__).parent.parent / 'data'
    ds = dups(p)
    assert not ds
