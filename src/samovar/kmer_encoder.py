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
import subprocess
import sys
import tempfile
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np

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


def _fmt(values: Iterable[float]) -> str:
    return "\t".join(f"{float(v):.6g}" for v in values)


def write_freq_table(
    r1: str,
    dest: str,
    k: int = DEFAULT_K,
    r2: Optional[str] = None,
    threads: int = 1,
) -> int:
    """Normalized k-mer frequencies via ``bin/samovar_kmer`` (one row per read)."""
    from samovar.kmer2 import ensure_kmer_binary

    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(ensure_kmer_binary()),
        "-i",
        r1,
        "-o",
        str(path),
        "-t",
        str(max(1, int(threads))),
        "-k",
        str(int(k)),
        "--normalize",
    ]
    if r2:
        cmd.extend(["-I", r2])
    subprocess.check_call(cmd)
    with path.open(encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def write_latent_table_stream(
    r1: str,
    r2: Optional[str],
    model_path: str,
    dest: str,
    device: str = "",
    threads: int = 1,
) -> int:
    """Encode C++ frequency rows in batches. A batch is dropped after it is written."""
    torch, model, dev, k, index, latent_dim = _load_encoder(model_path, device)
    names = [f"z{i}" for i in range(latent_dim)]
    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_in = len(index) if not isinstance(index, dict) else len(index)
    chunk = proba_chunk_rows(latent_dim, n_in)
    n = 0
    with tempfile.TemporaryDirectory(prefix="samovar-kmer-") as tmp:
        freq_path = Path(tmp) / "freq.tsv"
        write_freq_table(r1, str(freq_path), k=k, r2=r2, threads=threads)
        with freq_path.open(encoding="utf-8") as src, path.open("w", encoding="utf-8") as out:
            header = src.readline().rstrip("\n").split("\t")
            cols = header[1:]
            if isinstance(index, dict):
                places = [int(index[name]) for name in cols]
            else:
                places = list(range(len(cols)))
            out.write("seq\t" + "\t".join(names) + "\n")
            ids: List[str] = []
            vecs: List[np.ndarray] = []

            def _flush() -> None:
                nonlocal n
                if not ids:
                    return
                mat = np.zeros((len(vecs), n_in), dtype=np.float64)
                for row_i, vec in enumerate(vecs):
                    for src_i, dst_i in enumerate(places):
                        if 0 <= dst_i < n_in:
                            mat[row_i, dst_i] = vec[src_i]
                encoded = _encode_batch(torch, model, dev, mat)
                for read_id, row in zip(ids, encoded):
                    out.write(read_id + "\t" + _fmt(row) + "\n")
                    n += 1
                ids.clear()
                vecs.clear()

            for line in src:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                ids.append(parts[0])
                vecs.append(np.asarray(parts[1:], dtype=np.float64))
                if len(ids) >= chunk:
                    _flush()
            _flush()
    return n


def parse_output(file_path: str):
    from samovar.parse_annotators import read_custom_raw

    return read_custom_raw(file_path)


def extract_to_file(
    r1: str,
    r2: str,
    dest: str,
    db: str = "",
    k: int = DEFAULT_K,
    threads: int = 1,
) -> int:
    model = resolve_model_path(db)
    if model is not None:
        return write_latent_table_stream(
            r1, r2 or None, str(model), dest, threads=threads
        )
    return write_freq_table(r1, dest, k=k, r2=r2 or None, threads=threads)


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
    threads = int(args.threads)
    if threads < 1:
        parser.error("-t must be >= 1")
    extract_to_file(args.r1, args.r2, args.o, db=args.db, k=args.k, threads=threads)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
