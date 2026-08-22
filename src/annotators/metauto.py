#!/usr/bin/env python3

import argparse
import os
from collections import Counter
from itertools import product

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.preprocessing import LabelEncoder

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


from samovar.seqio import list_fasta_files, open_text, taxid_from_fasta_name


# ----------------------------
# FASTA / FASTQ readers
# ----------------------------
def read_fasta(path):
    with open_text(path) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    yield seq
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            yield seq


def read_fastq(path):
    with open_text(path) as f:
        while True:
            # Extract header to preserve the exact read ID required by parser_annotators
            h = f.readline().strip()
            if not h:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            # Return header (excluding '@') and sequence
            yield h[1:], seq


# ----------------------------
# k-mer encoding
# ----------------------------
def build_kmer_index(k):
    kmers = ["".join(p) for p in product("ACGT", repeat=k)]
    return {kmer: i for i, kmer in enumerate(kmers)}


def seq_to_vec(seq, k, kmer_index):
    vec = np.zeros(len(kmer_index))
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if kmer in kmer_index:
            vec[kmer_index[kmer]] += 1
    if vec.sum() > 0:
        vec /= vec.sum()
    return vec


# ----------------------------
# Model
# ----------------------------
class AutoEncoder(nn.Module):
    def __init__(self, input_dim, latent_dim, n_classes):
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
        x_rec = self.decoder(z)
        logits = self.classifier(z)
        return x_rec, logits, z


# ----------------------------
# TRAIN
# ----------------------------
def train(input_dir, out_model, k=4, epochs=10):

    kmer_index = build_kmer_index(k)

    X = []
    y = []

    print("[INFO] Loading genomes...")

    for path in list_fasta_files(input_dir, nucleotide=True, protein=False):
        taxid = taxid_from_fasta_name(path) or path.name
        for seq in read_fasta(str(path)):
            vec = seq_to_vec(seq.upper(), k, kmer_index)
            X.append(vec)
            y.append(taxid)

    X = np.array(X)
    le = LabelEncoder()
    y_enc = le.fit_transform(y)

    X = torch.tensor(X, dtype=torch.float32).to(DEVICE)
    y_enc = torch.tensor(y_enc).to(DEVICE)

    model = AutoEncoder(X.shape[1], latent_dim=64, n_classes=len(le.classes_)).to(
        DEVICE
    )

    opt = optim.Adam(model.parameters(), lr=1e-3)
    mse = nn.MSELoss()
    ce = nn.CrossEntropyLoss()

    print("[INFO] Training started...")

    for epoch in range(epochs):
        model.train()

        x_rec, logits, _ = model(X)

        loss_rec = mse(x_rec, X)
        loss_cls = ce(logits, y_enc)

        loss = loss_rec + loss_cls

        opt.zero_grad()
        loss.backward()
        opt.step()

        print(f"Epoch {epoch}: loss={loss.item():.4f}")

    torch.save(
        {
            "model": model.state_dict(),
            "k": k,
            "kmer_index": kmer_index,
            "label_encoder": le,
        },
        out_model,
    )

    print(f"[SUCCESS] Saved model: {out_model}")


# ----------------------------
# WORK
# ----------------------------
def work(fastq_dir, model_path):

    checkpoint = torch.load(model_path, map_location=DEVICE, weights_only=False)

    k = checkpoint["k"]
    kmer_index = checkpoint["kmer_index"]
    le = checkpoint["label_encoder"]

    model = AutoEncoder(len(kmer_index), 64, len(le.classes_)).to(DEVICE)
    model.load_state_dict(checkpoint["model"])
    model.eval()

    print("[INFO] Processing reads...")

    for fname in os.listdir(fastq_dir):
        if not (
            fname.endswith(".fastq") or fname.endswith(".fq") or fname.endswith(".gz")
        ):
            continue

        path = os.path.join(fastq_dir, fname)
        out_file = path + ".tsv"

        with open(out_file, "w") as out:
            # Write the exact header expected by parser_annotators.py
            out.write("seq\ttaxID\tconfidence\n")

            # Unpack header and seq from the updated read_fastq generator
            for header, seq in read_fastq(path):
                vec = seq_to_vec(seq.upper(), k, kmer_index)
                x = torch.tensor(vec, dtype=torch.float32).unsqueeze(0).to(DEVICE)

                with torch.no_grad():
                    _, logits, _ = model(x)
                    probs = torch.softmax(logits, dim=1).cpu().numpy()[0]

                idx = np.argmax(probs)
                taxid = le.inverse_transform([idx])[0]
                conf = probs[idx]

                # Write the actual read ID instead of a simple index
                out.write(f"{header}\t{taxid}\t{conf:.4f}\n")

        print(f"[SUCCESS] Saved predictions: {out_file}")


# ----------------------------
# CLI
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        prog="metauto", description="Metagenomic AutoEncoder Classifier"
    )

    sub = parser.add_subparsers(dest="cmd")

    t = sub.add_parser("train", help="Train the AutoEncoder model")
    t.add_argument("input_dir", help="Directory containing FASTA reference genomes")
    t.add_argument("out_model", help="Path to save the trained PyTorch model (.pt)")

    w = sub.add_parser("work", help="Classify reads using a trained model")
    w.add_argument("fastq_dir", help="Directory containing FASTQ files to classify")
    w.add_argument("model", help="Path to the trained PyTorch model")

    args = parser.parse_args()

    if args.cmd == "train":
        train(args.input_dir, args.out_model)
    elif args.cmd == "work":
        work(args.fastq_dir, args.model)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
