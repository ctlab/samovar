#!/usr/bin/env python3

import argparse
import os

import numpy as np
import torch

from samovar.kmer_encoder import (
    DEFAULT_K,
    DEFAULT_LATENT,
    make_autoencoder,
    seq_to_vec,
    train_encoder,
)
from samovar.seqio import open_text

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def read_fastq(path):
    with open_text(path) as f:
        while True:
            h = f.readline().strip()
            if not h:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            yield h[1:], seq


def train(input_dir, out_model, k=DEFAULT_K, epochs=10):
    train_encoder(input_dir, out_model, k=k, epochs=epochs, latent_dim=DEFAULT_LATENT)


def work(fastq_dir, model_path):
    checkpoint = torch.load(model_path, map_location=DEVICE, weights_only=False)
    k = checkpoint["k"]
    kmer_index = checkpoint["kmer_index"]
    le = checkpoint["label_encoder"]
    latent = int(checkpoint.get("latent_dim") or DEFAULT_LATENT)
    model = make_autoencoder(len(kmer_index), latent, len(le.classes_)).to(DEVICE)
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
            out.write("seq\ttaxID\tconfidence\n")
            for header, seq in read_fastq(path):
                vec = seq_to_vec(seq.upper(), k, kmer_index)
                x = torch.tensor(vec, dtype=torch.float32).unsqueeze(0).to(DEVICE)
                with torch.no_grad():
                    _, logits, _z = model(x)
                    probs = torch.softmax(logits, dim=1).cpu().numpy()[0]
                idx = np.argmax(probs)
                taxid = le.inverse_transform([idx])[0]
                conf = probs[idx]
                out.write(f"{header}\t{taxid}\t{conf:.4f}\n")
        print(f"[SUCCESS] Saved predictions: {out_file}")


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
