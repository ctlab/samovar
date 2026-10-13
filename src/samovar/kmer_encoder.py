#!/usr/bin/env python3
"""k-mer encoder Feature extractor (the encoding stage used inside metauto).

Without a trained model this writes normalized ACGT k-mer frequencies
(default k=4 → 256 columns). With ``-d model.pt`` it runs the autoencoder
and writes the latent vector (default 64 ``z*`` columns).

Train a model from reference FASTA (same algorithm as ``metauto train``):

    python -m samovar.kmer_encoder train -c GENOMES_DIR -o model.pt
"""

from __future__ import annotations

import argparse
import sys
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from samovar.kmer2 import iter_read_mates
from samovar.memory import proba_chunk_rows

DEFAULT_K = 4
DEFAULT_LATENT = 64
DEFAULT_EPOCHS = 10


def kmer_ids(k: int) -> List[str]:
    return ["".join(p) for p in product("ACGT", repeat=int(k))]


def build_kmer_index(k: int) -> Dict[str, int]:
    return {mer: i for i, mer in enumerate(kmer_ids(k))}


def seq_to_vec(seq: str, k: int, index: Optional[Dict[str, int]] = None) -> np.ndarray:
    """Normalized ACGT k-mer frequencies (metauto encoding)."""
    table = index or build_kmer_index(k)
    vec = np.zeros(len(table), dtype=np.float64)
    seq = (seq or "").upper()
    if k <= 0 or len(seq) < k:
        return vec
    for i in range(len(seq) - k + 1):
        mer = seq[i : i + k]
        idx = table.get(mer)
        if idx is not None:
            vec[idx] += 1
    total = vec.sum()
    if total > 0:
        vec /= total
    return vec


def _counts(seq: str, k: int, index: Dict[str, int]) -> np.ndarray:
    raw = np.zeros(len(index), dtype=np.float64)
    if k > 0 and seq and len(seq) >= k:
        for i in range(len(seq) - k + 1):
            idx = index.get(seq[i : i + k])
            if idx is not None:
                raw[idx] += 1
    return raw


def _pair_freq(left: str, right: str, k: int, index: Dict[str, int]) -> np.ndarray:
    """Normalized k-mer frequencies of R1 and R2, without a junction k-mer."""
    raw = _counts(left, k, index)
    if right:
        raw = raw + _counts(right, k, index)
    total = raw.sum()
    if total > 0:
        raw /= total
    return raw


def freqs_from_fastq(r1: str, r2: Optional[str] = None, k: int = DEFAULT_K) -> Dict[str, np.ndarray]:
    """In-memory frequencies. Sample-sized FASTQ should use :func:`write_freq_table`."""
    index = build_kmer_index(k)
    out: Dict[str, np.ndarray] = {}
    for read_id, left, right in iter_read_mates(r1, r2 or None):
        if not read_id:
            continue
        vec = _pair_freq(left, right, k, index)
        if read_id not in out:
            out[read_id] = vec
        else:
            merged = out[read_id] + vec
            total = merged.sum()
            out[read_id] = merged / total if total > 0 else merged
    return out


def resolve_model_path(db: str) -> Optional[Path]:
    raw = (db or "").strip()
    if not raw or raw in {".", "-"}:
        return None
    path = Path(raw).expanduser()
    if path.is_file():
        return path
    if path.is_dir():
        for name in ("kmer_encoder.pt", "model.pt", "metauto.pt"):
            cand = path / name
            if cand.is_file():
                return cand
        pts = sorted(path.glob("*.pt")) + sorted(path.glob("*.pth"))
        if pts:
            return pts[0]
    return None


def _torch_mods():
    import torch
    import torch.nn as nn
    import torch.optim as optim

    return torch, nn, optim


def make_autoencoder(input_dim: int, latent_dim: int, n_classes: int):
    _torch, nn, _optim = _torch_mods()

    class AutoEncoder(nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = nn.Sequential(
                nn.Linear(input_dim, 512), nn.ReLU(), nn.Linear(512, latent_dim)
            )
            self.decoder = nn.Sequential(
                nn.Linear(latent_dim, 512), nn.ReLU(), nn.Linear(512, input_dim)
            )
            self.classifier = nn.Sequential(
                nn.Linear(latent_dim, 128), nn.ReLU(), nn.Linear(128, n_classes)
            )

        def forward(self, x):
            z = self.encoder(x)
            return self.decoder(z), self.classifier(z), z

    return AutoEncoder()


def train_encoder(
    input_dir: str | Sequence[str],
    dest: str,
    k: int = DEFAULT_K,
    epochs: int = DEFAULT_EPOCHS,
    latent_dim: int = DEFAULT_LATENT,
    device: str = "",
) -> Path:
    """Train reconstruction + taxid classifier on k-mer frequencies of FASTA genomes."""
    torch, nn, optim = _torch_mods()
    from sklearn.preprocessing import LabelEncoder

    from samovar.seqio import list_fasta_files, open_text, taxid_from_fasta_name

    dirs = [input_dir] if isinstance(input_dir, (str, Path)) else list(input_dir)
    index = build_kmer_index(k)
    X: List[np.ndarray] = []
    y: List[str] = []
    for folder in dirs:
        for path in list_fasta_files(folder, nucleotide=True, protein=False):
            taxid = str(taxid_from_fasta_name(path) or path.stem)
            seq = ""
            with open_text(str(path)) as handle:
                for line in handle:
                    if line.startswith(">"):
                        if seq:
                            X.append(seq_to_vec(seq, k, index))
                            y.append(taxid)
                            seq = ""
                    else:
                        seq += line.strip()
                if seq:
                    X.append(seq_to_vec(seq, k, index))
                    y.append(taxid)
    if not X:
        raise ValueError(f"no nucleotide FASTA under {dirs}")
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    dev = device or ("cuda" if torch.cuda.is_available() else "cpu")
    xt = torch.tensor(np.array(X), dtype=torch.float32).to(dev)
    yt = torch.tensor(y_enc).to(dev)
    model = make_autoencoder(xt.shape[1], latent_dim, len(le.classes_)).to(dev)
    opt = optim.Adam(model.parameters(), lr=1e-3)
    mse = nn.MSELoss()
    ce = nn.CrossEntropyLoss()
    for epoch in range(int(epochs)):
        model.train()
        rec, logits, _z = model(xt)
        loss = mse(rec, xt) + ce(logits, yt)
        opt.zero_grad()
        loss.backward()
        opt.step()
        print(f"Epoch {epoch}: loss={float(loss.item()):.4f}", file=sys.stderr)
    out = Path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "model": model.state_dict(),
            "k": int(k),
            "kmer_index": index,
            "label_encoder": le,
            "latent_dim": int(latent_dim),
        },
        out,
    )
    return out


def _load_encoder(model_path: str, device: str = ""):
    torch, _nn, _optim = _torch_mods()
    ckpt = torch.load(model_path, map_location="cpu", weights_only=False)
    k = int(ckpt["k"])
    index = ckpt["kmer_index"]
    le = ckpt["label_encoder"]
    latent_dim = int(ckpt.get("latent_dim") or 64)
    dev = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model = make_autoencoder(len(index), latent_dim, len(le.classes_)).to(dev)
    model.load_state_dict(ckpt["model"])
    model.eval()
    return torch, model, dev, k, index, latent_dim


def _encode_batch(torch, model, dev, mat: np.ndarray) -> np.ndarray:
    with torch.no_grad():
        xt = torch.tensor(mat, dtype=torch.float32).to(dev)
        _rec, _logits, z = model(xt)
        return z.detach().cpu().numpy()


def latents_from_fastq(
    r1: str,
    r2: Optional[str],
    model_path: str,
    device: str = "",
) -> Tuple[List[str], np.ndarray]:
    """Encode in memory-sized batches. The returned matrix is the caller's."""
    torch, model, dev, k, index, latent_dim = _load_encoder(model_path, device)
    ids: List[str] = []
    blocks: List[np.ndarray] = []
    batch_ids: List[str] = []
    batch_vecs: List[np.ndarray] = []
    chunk = proba_chunk_rows(latent_dim, len(index))

    def _flush() -> None:
        if not batch_ids:
            return
        blocks.append(_encode_batch(torch, model, dev, np.stack(batch_vecs)))
        ids.extend(batch_ids)
        batch_ids.clear()
        batch_vecs.clear()

    for read_id, left, right in iter_read_mates(r1, r2 or None):
        if not read_id:
            continue
        batch_ids.append(read_id)
        batch_vecs.append(_pair_freq(left, right, k, index))
        if len(batch_ids) >= chunk:
            _flush()
    _flush()
    if not ids:
        return [], np.zeros((0, latent_dim))
    return ids, np.vstack(blocks)


def write_latent_table_stream(
    r1: str,
    r2: Optional[str],
    model_path: str,
    dest: str,
    device: str = "",
) -> int:
    """Write latent rows in batches. Each batch is released after it is written."""
    torch, model, dev, k, index, latent_dim = _load_encoder(model_path, device)
    names = [f"z{i}" for i in range(latent_dim)]
    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    chunk = proba_chunk_rows(latent_dim, len(index))
    n = 0
    batch_ids: List[str] = []
    batch_vecs: List[np.ndarray] = []
    with path.open("w", encoding="utf-8") as out:
        out.write("seq\t" + "\t".join(names) + "\n")

        def _flush() -> None:
            nonlocal n
            if not batch_ids:
                return
            encoded = _encode_batch(torch, model, dev, np.stack(batch_vecs))
            for read_id, row in zip(batch_ids, encoded):
                out.write(read_id + "\t" + _fmt(row) + "\n")
                n += 1
            batch_ids.clear()
            batch_vecs.clear()

        for read_id, left, right in iter_read_mates(r1, r2 or None):
            if not read_id:
                continue
            batch_ids.append(read_id)
            batch_vecs.append(_pair_freq(left, right, k, index))
            if len(batch_ids) >= chunk:
                _flush()
        _flush()
    return n


def _fmt(values: Iterable[float]) -> str:
    return "\t".join(f"{float(v):.6g}" for v in values)


def write_freq_table(r1: str, dest: str, k: int = DEFAULT_K, r2: Optional[str] = None) -> int:
    """Write one frequency row per read. Counters for a read are dropped after the row."""
    names = kmer_ids(k)
    index = build_kmer_index(k)
    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8") as out:
        out.write("seq\t" + "\t".join(names) + "\n")
        for read_id, left, right in iter_read_mates(r1, r2 or None):
            if not read_id:
                continue
            out.write(read_id + "\t" + _fmt(_pair_freq(left, right, k, index)) + "\n")
            n += 1
    return n


def write_latent_table(ids: Sequence[str], z: np.ndarray, dest: str) -> int:
    n = 0 if z.size == 0 else z.shape[1]
    names = [f"z{i}" for i in range(n)]
    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        out.write("seq\t" + "\t".join(names) + "\n")
        for i, read_id in enumerate(ids):
            out.write(read_id + "\t" + _fmt(z[i]) + "\n")
    return len(ids)


def parse_output(file_path: str):
    from samovar.parse_annotators import read_custom_raw

    return read_custom_raw(file_path)


def extract_to_file(r1: str, r2: str, dest: str, db: str = "", k: int = DEFAULT_K) -> int:
    model = resolve_model_path(db)
    if model is not None:
        return write_latent_table_stream(r1, r2 or None, str(model), dest)
    return write_freq_table(r1, dest, k=k, r2=r2 or None)


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in {"train", "encode", "freqs"}:
        cmd = argv[0]
        rest = argv[1:]
    else:
        cmd = "extract"
        rest = argv

    if cmd == "train":
        parser = argparse.ArgumentParser(prog="python -m samovar.kmer_encoder train")
        parser.add_argument("-c", "--input-dir", dest="input_dir", required=True)
        parser.add_argument("-o", dest="o", required=True)
        parser.add_argument("-k", dest="k", type=int, default=DEFAULT_K)
        parser.add_argument("--epochs", type=int, default=DEFAULT_EPOCHS)
        parser.add_argument("--latent", type=int, default=DEFAULT_LATENT)
        args = parser.parse_args(rest)
        train_encoder(args.input_dir, args.o, k=args.k, epochs=args.epochs, latent_dim=args.latent)
        return 0

    parser = argparse.ArgumentParser(
        description="k-mer encoder Feature table (frequencies or autoencoder latents)."
    )
    parser.add_argument("-i", "-1", dest="r1", required=True)
    parser.add_argument("-I", "-2", dest="r2", default="")
    parser.add_argument("-d", dest="db", default="")
    parser.add_argument("-o", "--output", dest="o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    parser.add_argument("-k", dest="k", type=int, default=DEFAULT_K)
    args = parser.parse_args(rest)
    extract_to_file(args.r1, args.r2, args.o, db=args.db, k=args.k)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
