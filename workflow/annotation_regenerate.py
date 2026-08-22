"""CLI for metagenome annotation table regeneration (Python path)."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from samovar.regenerate import normalize_regeneration_mode
from samovar.table2iss import samovar_annotation_regenerate


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Regenerate per-annotator abundance tables from annotation CSVs."
    )
    parser.add_argument("--annotation_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--config", help="YAML config (samovaR-style keys)")
    parser.add_argument(
        "--regeneration_mode",
        default="direct",
        choices=["direct", "preserve", "glm", "bootstrap", "vae", "samovar"],
    )
    parser.add_argument(
        "--N",
        type=int,
        default=None,
        help="Synthetic sample count (default: same as input samples)",
    )
    parser.add_argument("--N_reads", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--rescale_abundance",
        action="store_true",
        help="Force synthetic column totals to N_reads (off by default).",
    )
    args = parser.parse_args()

    cfg_path = args.config
    if cfg_path:
        with open(cfg_path, encoding="utf-8") as handle:
            cfg = yaml.safe_load(handle) or {}
    else:
        cfg = {}
    cfg["regeneration_mode"] = normalize_regeneration_mode(
        args.regeneration_mode or cfg.get("regeneration_mode", "direct")
    )
    if args.N is not None:
        cfg["N"] = args.N
    cfg["N_reads"] = args.N_reads
    cfg["seed"] = args.seed
    cfg["rescale_abundance"] = bool(args.rescale_abundance)
    cfg["output_dir"] = args.output_dir

    tmp = Path(args.output_dir) / ".regeneration_config.yaml"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    with open(tmp, "w", encoding="utf-8") as handle:
        yaml.dump(cfg, handle)

    samovar_annotation_regenerate(
        annotation_dir=args.annotation_dir,
        config_samovar=str(tmp),
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
