"""Stable MAG ids from tetranucleotide homology, and their use by RandomForest.

MAG assembly is outside the Kraken2/Kaiju install path, so this module is optional.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytestmark = pytest.mark.optional

from samovar.MAG.homology import assign_ids, stabilize_mags
from samovar.MAG.mapping import map_reads
from samovar.MAG.tetra import DEFAULT_THRESHOLD, kmer_counts, profile_rms_z, read_fasta_sequence
from samovar.reprofiling import train_models

ROOT = Path(__file__).resolve().parents[1]
TOY = [
    ROOT / "data/test_genomes/meta/562.fna",
    ROOT / "data/test_genomes/meta/4932.fna",
    ROOT / "data/test_genomes/host/9606.fna",
    ROOT / "data/test_genomes/meta/2886930.fna",
]
# Optional local phage example. Not in git; GCF_000819615 matches bundled 2886930.
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
    # The yeast-only reassembly must keep yeast's id, not become M0001 by recounting.
    yeast_only = tmp_path / "yeast_only"
    _write(yeast_only, "other_name.fa", _mutate(yeast, seed=4))
    yeast_ids = _ids_for(yeast_only, catalog, tmp_path / "stable_yeast")
    assert yeast_ids["other_name"] == first["bin3"]
    assert yeast_ids["other_name"] != first["bin17"]
    second = _ids_for(regen, catalog, tmp_path / "stable_regen")
    assert second["contig_zz"] == first["bin17"]
    assert second["contig_aa"] == first["bin3"]


def test_toy_and_phage_profiles_separate_genomes():
    sequences = {path.name: _slice(path) for path in TOY}
    phage_files = sorted(PHAGE.glob("*.fa.gz")) if PHAGE.is_dir() else []
    # GCF_000819615 is the same short sequence as bundled 2886930.
    for path in phage_files:
        if path.name.startswith("GCF_000819615"):
            continue
        sequences[path.name] = _slice(path)
        if len(sequences) >= 6:
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


def _fastq(path: Path, records) -> None:
    lines = []
    for name, sequence in records:
        lines.extend([f"@{name}", sequence, "+", "I" * len(sequence)])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_mapping_onto_initial_mags_on_toy_and_phage(tmp_path):
    """Reads from a regenerated copy keep the initial MAG id, and RandomForest uses it."""
    sequences = {path.stem: _slice(path) for path in TOY}
    phage_files = sorted(PHAGE.glob("*.fa.gz")) if PHAGE.is_dir() else []
    for path in phage_files:
        if path.name.startswith("GCF_000819615"):
            continue
        sequences[path.stem] = _slice(path)
        break
    assert len(sequences) >= 2
    assert phage_files or len(sequences) >= 2
    initial = tmp_path / "initial_mags"
    catalog = tmp_path / "MAG" / "catalog.tsv"
    names = list(sequences)
    for i, name in enumerate(names):
        _write(initial, f"bin_{i}.fa", sequences[name])
    stable = _ids_for(initial, catalog, tmp_path / "stable")
    by_mag = {sequences[name]: stable[f"bin_{i}"] for i, name in enumerate(names)}
    records = []
    truth = []
    taxa = {names[0]: 562, names[1]: 4932}
    for i, name in enumerate(names):
        read = sequences[name][100:180]
        records.append((f"r{i}", read))
        truth.append(taxa.get(name, 10000 + i))
    fastq = tmp_path / "r.fastq"
    _fastq(fastq, records)
    placed = {read_id: mag for read_id, _seq, mag in map_reads(str(fastq), "", tmp_path / "stable", fallback=False)}
    for i, name in enumerate(names):
        assert placed[f"r{i}"] == by_mag[sequences[name]]

    rng = np.random.default_rng(1)
    n = 60
    owners = [names[i % len(names)] for i in range(n)]
    true = np.array([taxa.get(owner, 10000 + names.index(owner)) for owner in owners])
    pool = sorted(set(int(v) for v in true))
    noisy = true.copy()
    for i in range(n):
        if rng.random() < 0.5:
            choices = [tax for tax in pool if tax != int(true[i])]
            noisy[i] = choices[int(rng.integers(0, len(choices)))]
    mag = [by_mag[sequences[owner]] for owner in owners]
    frame = pd.DataFrame(
        {
            "seq": [f"q{i}" for i in range(n)],
            "taxid_dummy": noisy,
            "feat_assembly_0_MAG_ID": mag,
            "length": 80,
            "true": true,
        }
    )
    _model, _models, metrics, features = train_models(
        frame, test_size=0.25, random_state=0, methods="random_forest"
    )
    assert "feat_assembly_0_MAG_ID" in features
    bare = frame.drop(columns=["feat_assembly_0_MAG_ID"])
    _model, _models, bare_metrics, _features = train_models(
        bare, test_size=0.25, random_state=0, methods="random_forest"
    )
    assert metrics["RandomForest"] > bare_metrics["RandomForest"] + 0.15


def test_assign_ids_does_not_use_filename(tmp_path):
    left = tmp_path / "a"
    right = tmp_path / "b"
    seq = "ACGT" * 500
    _write(left, "random_9.fa", seq)
    _write(right, "random_1.fa", seq)
    rows_left = assign_ids(list(left.glob("*.fa")), [])
    rows_right = assign_ids(list(right.glob("*.fa")), [])
    assert rows_left[0]["mag_id"] == rows_right[0]["mag_id"]
