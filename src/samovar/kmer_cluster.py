#!/usr/bin/env python3
"""Train and apply a shared k-mer clustering feature model."""
from __future__ import annotations

import argparse
import itertools
import os
import random
import tempfile
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

from samovar.kmer2 import iter_read_mates
from samovar.kmer_encoder import kmer_ids, write_freq_table


def _vectors(records, k, threads, scratch, reverse_complemented=True):
    """Use the existing native encoder on bounded FASTQ batches."""
    if not records:
        return [], np.empty((0, len(kmer_ids(k, reverse_complemented))), dtype=np.float32)
    r1, r2, table = (Path(scratch) / name for name in ('r1.fastq', 'r2.fastq', 'freq.tsv'))
    # Unique synthetic ids prevent collisions during native pair merging.
    with r1.open('w') as first, r2.open('w') as second:
        for i, (_, a, b) in enumerate(records):
            for handle, seq in ((first, a), (second, b)):
                handle.write(f'@{i}\n{seq}\n+\n{"I" * len(seq)}\n')
    write_freq_table(str(r1), str(table), k=k, r2=str(r2), threads=threads,
                     reverse_complemented=reverse_complemented)
    frame = pd.read_csv(table, sep='\t', dtype={'seq': str})
    frame['seq'] = pd.to_numeric(frame['seq'], errors='raise').astype(int)
    frame = frame.set_index('seq').reindex(range(len(records)))
    values = frame[kmer_ids(k, reverse_complemented)].fillna(0).to_numpy(dtype=np.float32)
    return [row[0] for row in records], values


def _batch_size(k, memory_mb=256, reverse_complemented=True):
    width = len(kmer_ids(k, reverse_complemented))
    return max(1, min(4096, memory_mb * 1024 * 1024 // max(1, width * 16)))


def _reservoir(inputs, limit, seed, mode):
    rng = random.Random(seed)
    samples = []
    if mode == 'balanced':
        if limit < len(inputs):
            raise ValueError('Training limit is smaller than the number of initial samples')
        # Independently bounded reservoirs, then distribute unused slots fairly.
        limit_each = max(1, limit // len(inputs))
        for r1, r2 in inputs:
            reservoir = []
            for n, record in enumerate(iter_read_mates(r1, r2), 1):
                if len(reservoir) < limit_each:
                    reservoir.append(record)
                else:
                    pos = rng.randrange(n)
                    if pos < limit_each:
                        reservoir[pos] = record
            samples.extend(reservoir)
        return samples[:limit]
    for n, record in enumerate(itertools.chain.from_iterable(iter_read_mates(a, b) for a, b in inputs), 1):
        if len(samples) < limit:
            samples.append(record)
        else:
            pos = rng.randrange(n)
            if pos < limit:
                samples[pos] = record
    return samples


def train(inputs, dest, *, k=4, min_clusters=1, max_clusters=20, max_reads=50000,
          sampling='proportional', seed=0, threads=1, memory_mb=256, reverse_complemented=True):
    if not 1 <= k <= 8 or not 1 <= min_clusters <= max_clusters or max_reads < 1 or not inputs:
        raise ValueError('Invalid k, cluster range, training size, or empty input cohort')
    if sampling not in ('proportional', 'balanced'):
        raise ValueError('sampling must be proportional or balanced')
    # Bound the training matrix and the temporary frequency batch independently.
    limit = min(max_reads, max(1, memory_mb * 1024 * 1024 // (len(kmer_ids(k, reverse_complemented)) * 12)))
    records = _reservoir(inputs, limit, seed, sampling)
    batch = _batch_size(k, memory_mb, reverse_complemented)
    matrices = []
    with tempfile.TemporaryDirectory(prefix='samovar-kmer-cluster-') as tmp:
        for offset in range(0, len(records), batch):
            _, values = _vectors(records[offset:offset + batch], k, threads, tmp, reverse_complemented)
            matrices.append(values[values.sum(axis=1) > 0])
    x = np.concatenate(matrices) if matrices else np.empty((0, len(kmer_ids(k, reverse_complemented))), dtype=np.float32)
    if not len(x):
        raise ValueError('Cannot train kmer-cluster: no reads contain a valid k-mer')
    unique = len(np.unique(x, axis=0))
    high = min(max_clusters, unique, len(x))
    low = min(min_clusters, high)
    models = []
    for n in range(low, high + 1):
        model = KMeans(n_clusters=n, random_state=seed, n_init=10).fit(x)
        models.append(model)
    inertias = np.array([m.inertia_ for m in models], dtype=float)
    selected = 0
    if len(models) > 2 and inertias[0] > inertias[-1]:
        xs = np.linspace(0, 1, len(models))
        ys = (inertias - inertias[-1]) / (inertias[0] - inertias[-1])
        distances = (1 - xs) - ys
        selected = int(np.argmax(distances))
        if distances[selected] < 0.02:
            warnings.warn("Weak k-mer clustering elbow; using the best available candidate")
    elif len(models) <= 2:
        warnings.warn("Too few candidate cluster counts for an elbow; using the first candidate")
    fingerprints = [(Path(path).name, Path(path).stat().st_size, Path(path).stat().st_mtime_ns)
                    for pair in inputs for path in pair if path]
    payload = {'version': 2, 'model': models[selected], 'k': k,
               'reverse_complemented': reverse_complemented, 'kmers': kmer_ids(k, reverse_complemented),
               'inputs': fingerprints, 'options': {'k': k, 'min_clusters': min_clusters,
                   'max_clusters': max_clusters, 'max_reads': max_reads, 'sampling': sampling,
                   'seed': seed, 'memory_mb': memory_mb, 'reverse_complemented': reverse_complemented},
               'inertias': inertias.tolist(), 'cluster_range': [low, high],
               'sampling': sampling, 'seed': seed, 'training_reads': len(x),
               'selected_clusters': models[selected].n_clusters}
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=dest.parent, prefix='.kmer-cluster-', delete=False) as tmp:
        name = tmp.name
    try:
        joblib.dump(payload, name)
        os.replace(name, dest)
    finally:
        Path(name).unlink(missing_ok=True)
    return payload


def predict(r1, r2, model_path, dest, *, threads=1, memory_mb=256):
    payload = joblib.load(model_path)
    if payload.get('version') != 2 or payload.get('kmers') != kmer_ids(payload['k'], payload['reverse_complemented']):
        raise ValueError('Incompatible kmer-cluster model')
    k = payload['k']
    batch = _batch_size(k, memory_mb, payload["reverse_complemented"])
    records = iter_read_mates(r1, r2 or None)
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='samovar-kmer-cluster-') as tmp, open(dest, 'w') as out:
        out.write('seq\tfeat_kmer_cluster\n')
        while chunk := list(itertools.islice(records, batch)):
            ids, x = _vectors(chunk, k, threads, tmp, payload["reverse_complemented"])
            valid = x.sum(axis=1) > 0
            labels = np.full(len(x), -1, dtype=int)
            if valid.any():
                labels[valid] = payload['model'].predict(x[valid])
            for name, label in zip(ids, labels):
                out.write(f'{name}\t{label if label >= 0 else ""}\n')


def parse_output(path):
    from samovar.parse_annotators import read_custom_raw
    return read_custom_raw(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', required=True)
    parser.add_argument('-I', default='')
    parser.add_argument('-d', required=True, help='Saved k-mer cluster model')
    parser.add_argument('-o', required=True)
    parser.add_argument('-t', type=int, default=1)
    parser.add_argument('--memory-mb', type=int, default=256)
    args = parser.parse_args(argv)
    predict(args.i, args.I, args.d, args.o, threads=args.t, memory_mb=args.memory_mb)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())


def training_options(extra):
    """Parse the feature's training flags from an annotator's extra options."""
    import shlex
    parser = argparse.ArgumentParser(prog='kmer-cluster options')
    parser.add_argument('-k', type=int, default=4)
    parser.add_argument('--min-clusters', type=int, default=1)
    parser.add_argument('--max-clusters', type=int, default=20)
    parser.add_argument('--max-reads', type=int, default=50000)
    parser.add_argument('--sampling', choices=('balanced', 'proportional'), default='proportional')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--memory-mb', type=int, default=256)
    parser.add_argument("--no-reverse-complement", dest="reverse_complemented", action="store_false", default=True)
    return vars(parser.parse_args(shlex.split(extra or '')))


def validate_model(path, extra, inputs=None):
    payload = joblib.load(path)
    if payload.get("version") != 2 or payload.get("options") != training_options(extra):
        raise ValueError("Incompatible kmer-cluster model settings; rebuild the initial annotations and dependent results")
    if inputs is not None:
        fingerprints = [(Path(p).name, Path(p).stat().st_size, Path(p).stat().st_mtime_ns)
                        for pair in inputs for p in pair if p]
        if payload.get("inputs") != fingerprints:
            raise ValueError("kmer-cluster initial inputs changed; rebuild the initial annotations and dependent results")
    return payload


def cluster_levels_for_run(root, columns):
    """Return every fitted cluster ID for categorical model preprocessing."""
    import yaml
    config = Path(root) / ".log" / "configs" / "config_init.yaml"
    if not config.is_file():
        return {}
    runs = (yaml.safe_load(config.read_text()) or {}).get("run_config") or []
    sizes = []
    for run in runs:
        if str(run.get("type", "")).replace("-", "_") != "kmer_cluster":
            continue
        db = str(run.get("db_path") or "")
        path = Path(db) if db not in ("", ".") else Path(root) / "kmer_cluster_models" / run["run_name"] / "model.joblib"
        payload = joblib.load(path)
        if payload.get("version") != 2:
            raise ValueError(f"Incompatible kmer-cluster model: {path}")
        sizes.append(payload["model"].n_clusters)
    # Multiple runs may be ingested in a different order from the config.
    # The union only adds all-zero columns when their sizes differ.
    if not sizes:
        return {}
    return {col: list(range(max(sizes))) for col in columns
            if str(col).startswith("feat_kmer_cluster_") and str(col).endswith("_kmer_cluster")}
