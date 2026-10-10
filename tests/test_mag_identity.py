"""Stable MAG ids from tetranucleotide homology, and their use by RandomForest."""

from pathlib import Path

import numpy as np
import pandas as pd

from samovar.MAG.homology import assign_ids, stabilize_mags
from samovar.MAG.tetra import DEFAULT_THRESHOLD, kmer_counts, profile_rms_z, read_fasta_sequence
from samovar.reprofiling import train_models

ROOT = Path(__file__).resolve().parents[1]
TOY = [
    ROOT / "data/test_genomes/meta/562.fna",
    ROOT / "data/test_genomes/meta/4932.fna",
]
PHAGE = (
    ROOT
    / "examples_outdir/phage/database1/.database/kaiju_src/.genomes/processed"
)


def _slice(path: Path, n: int = 12000) -> str:
    return read_fasta_sequence(path)[:n]


def _write(dir_path: Path, name: str, sequence: str) -> None:
    dir_path.mkdir(parents=True, exist_ok=True)
    wrapped = "\n".join(sequence[i : i + 80] for i in range(0, len(sequence), 80))
    (dir_path / name).write_text(f">{name}\n{wrapped}\n", encoding="utf-8")


def _mutate(sequence: str, rate: float = 0.01, seed: int = 1) -> str:
    rng = np.random.default_rng(seed)
    bases = np.array(list("ACGT"))
    out = list(sequence)
    for i, base in enumerate(out):
        if rng.random() < rate:
            out[i] = str(rng.choice(bases[bases != base]))
    return "".join(out)


def _ids_for(directory: Path, catalog: Path, fasta_out: Path) -> dict:
    rows = stabilize_mags(directory, catalog, fasta_out)
    return {row["query"]: row["mag_id"] for row in rows}


def test_regenerated_bins_reuse_initial_ids(tmp_path):
    initial = tmp_path / "initial"
    regen = tmp_path / "regen"
    catalog = tmp_path / "MAG" / "catalog.tsv"
    ecoli = _slice(TOY[0])
    yeast = _slice(TOY[1])
    _write(initial, "bin17.fa", ecoli)
    _write(initial, "bin3.fa", yeast)
    first = _ids_for(initial, catalog, tmp_path / "stable_initial")
    assert first["bin17"] != first["bin3"]
    _write(regen, "contig_zz.fa", _mutate(ecoli, seed=2))
    _write(regen, "contig_aa.fa", _mutate(yeast, seed=3))
    second = _ids_for(regen, catalog, tmp_path / "stable_regen")
    assert second["contig_zz"] == first["bin17"]
    assert second["contig_aa"] == first["bin3"]


def test_toy_and_phage_profiles_separate_genomes():
    sequences = {path.name: _slice(path) for path in TOY}
    phage_files = sorted(PHAGE.glob("*.fa.gz")) if PHAGE.is_dir() else []
    # GCF_000819615 is the same short sequence as toy 2886930; keep a distinct phage.
    for path in phage_files:
        if path.name.startswith("GCF_000819615"):
            continue
        sequences[path.name] = _slice(path)
        if len(sequences) >= 4:
            break
    assert len(sequences) >= 3
    names = list(sequences)
    for i, left in enumerate(names):
        mutated = _mutate(sequences[left], seed=10 + i)
        counts_a, n_a = kmer_counts(sequences[left])
        counts_b, n_b = kmer_counts(mutated)
        assert profile_rms_z(counts_a, n_a, counts_b, n_b) < DEFAULT_THRESHOLD
        for right in names[i + 1 :]:
            counts_c, n_c = kmer_counts(sequences[right])
            assert profile_rms_z(counts_a, n_a, counts_c, n_c) > DEFAULT_THRESHOLD


def test_random_forest_gains_from_stable_mag_id():
    rng = np.random.default_rng(0)
    n = 80
    true = np.array([562] * (n // 2) + [9606] * (n // 2))
    noisy = true.copy()
    flip = rng.random(n) < 0.45
    noisy[flip] = np.where(noisy[flip] == 562, 9606, 562)
    stable = np.where(true == 562, "M0001", "M0002")
    shuffled = stable.copy()
    rng.shuffle(shuffled)

    def _score(mag):
        frame = pd.DataFrame(
            {
                "seq": [f"r{i}" for i in range(n)],
                "taxid_dummy": noisy,
                "feat_assembly_0_MAG_ID": mag,
                "length": 100,
                "true": true,
            }
        )
        _model, _models, metrics, features = train_models(
            frame, test_size=0.25, random_state=0, methods="random_forest"
        )
        assert "feat_assembly_0_MAG_ID" in features
        return metrics["RandomForest"]

    with_mag = _score(stable)
    without = _score(["NA"] * n)
    scrambled = _score(shuffled)
    assert with_mag > without + 0.15
    assert with_mag > scrambled + 0.15


def test_assign_ids_does_not_use_filename(tmp_path):
    left = tmp_path / "a"
    right = tmp_path / "b"
    seq = "ACGT" * 500
    _write(left, "random_9.fa", seq)
    _write(right, "random_1.fa", seq)
    rows_left = assign_ids(list(left.glob("*.fa")), [])
    rows_right = assign_ids(list(right.glob("*.fa")), [])
    assert rows_left[0]["mag_id"] == rows_right[0]["mag_id"]
