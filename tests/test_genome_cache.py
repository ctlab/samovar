"""Bundled test_genomes are never an NCBI library. Index/reuse lives in test_reindex.py."""

from __future__ import annotations

from samovar.genome_cache import is_bundled_test_genomes_path
from samovar.paths import test_genomes_dir as bundled_test_genomes_dir


def test_refuses_bundled_test_genomes_as_library():
    tg = bundled_test_genomes_dir()
    assert is_bundled_test_genomes_path(tg)
    assert is_bundled_test_genomes_path(tg / "meta")
    assert is_bundled_test_genomes_path(tg / "host" / "9606.fna")
    assert not is_bundled_test_genomes_path(tg.parent)
