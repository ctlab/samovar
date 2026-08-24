#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

trap 'rm -rf "$TMP_DIR"' EXIT
# Default threads
THREADS=4

# Parse command line arguments
while getopts "i:I:d:o:p:t:" opt; do
  case $opt in
    i) R1="$OPTARG" ;;
    I) R2="$OPTARG" ;;
    d) DB="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    p) TOOL="$OPTARG" ;; # Used as tool selector (e.g., "metauto", "centrifuge")
    t) THREADS="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Validate required arguments
if [ -z "$R1" ] || [ -z "$R2" ] || [ -z "$OUT" ] || [ -z "$TOOL" ]; then
  echo "[ERROR] Missing required arguments: -i, -I, -o, or -p"
  echo "Usage: ./custom.sh -i R1 -I R2 -d DB -o OUT -p TOOL_NAME [-t THREADS]"
  exit 1
fi

# Create a temporary working directory for this specific run
TMP_DIR="$(dirname "$OUT")/tmp_$(basename "$OUT" .out)"
mkdir -p "$TMP_DIR"

# Check if input FASTQ is virtually empty to prevent downstream tool crashes
if [[ "$R1" == *.gz ]]; then
  LINES=$(gzip -dc "$R1" | wc -l)
else
  LINES=$(wc -l < "$R1")
fi
if [ "$LINES" -lt 4 ]; then
    echo "[WARNING] Input fastq is empty or malformed. Generating empty output."
    touch "$OUT"
    exit 0
fi

echo "[INFO] Running custom wrapper for tool: ${TOOL}"

# ==========================================
# TOOL ROUTER
# ==========================================
case "$TOOL" in

  # ----------------------------------------
  # Dummy / constant-taxID classifier (pytest + custom-class smoke test)
  # ----------------------------------------
  "dummy"|"dummy9606"|"constant9606"|"constant"|"random")
    echo "[INFO] Starting constant-taxID dummy classifier (9606)..."
    python "$SCRIPT_DIR/constant9606.py" -i "$R1" -I "$R2" -o "$OUT" --taxid 9606
    ;;

  # ----------------------------------------
  # 1. Deep Learning AutoEncoder (metauto)
  # ----------------------------------------
  "metauto")
    echo "[INFO] Starting metauto DL classification..."
    # metauto expects a directory of FASTQ files. We symlink them.
    FQ_DIR="$TMP_DIR/fastq"
    mkdir -p "$FQ_DIR"
    if [[ "$R1" == *.gz ]]; then
      ln -sf "$(realpath "$R1")" "$FQ_DIR/R1.fastq.gz"
      ln -sf "$(realpath "$R2")" "$FQ_DIR/R2.fastq.gz"
    else
      ln -sf "$(realpath "$R1")" "$FQ_DIR/R1.fastq"
      ln -sf "$(realpath "$R2")" "$FQ_DIR/R2.fastq"
    fi

    # Run the python script
    python "$SCRIPT_DIR/metauto.py" work "$FQ_DIR" "$DB"

    # Concat output TSVs, remove headers, and take only seq and taxID (columns 1 and 2)
    cat "$FQ_DIR"/*.tsv | grep -v "^seq" | awk -F'\t' '{print $1 "\t" $2}' > "$OUT"
    ;;

  # ----------------------------------------
  # 2. Centrifuge Classifier
  # ----------------------------------------
  "centrifuge")
    echo "[INFO] Starting Centrifuge classification..."
    # -x wants an index prefix. Accept a directory of *.1.cf, a .1.cf path, or a prefix.
    INDEX="$DB"
    if [ -d "$DB" ]; then
      cf1=$(ls "$DB"/*.1.cf 2>/dev/null | head -n 1 || true)
      if [ -n "$cf1" ]; then
        INDEX="${cf1%.1.cf}"
      fi
    elif [[ "$DB" == *.1.cf ]]; then
      INDEX="${DB%.1.cf}"
    fi
    echo "[INFO] Centrifuge index prefix: ${INDEX}"
    # --report-file: centrifuge otherwise writes centrifuge_report.tsv into CWD
    # (often the repo root when snakemake cds there). Keep all litter in TMP_DIR.
    centrifuge -x "$INDEX" -1 "$R1" -2 "$R2" -p "$THREADS" \
      -S "$TMP_DIR/centrifuge.out" \
      --report-file "$TMP_DIR/centrifuge_report.tsv"
    
    # Parse output: skip header (NR>1), only print classified reads (taxID != 0)
    awk -F'\t' '{if ($3 != 0 && NR>1) print $1 "\t" $3}' "$TMP_DIR/centrifuge.out" > "$OUT"
    ;;

  # ----------------------------------------
  # 3. Assembly + Mapping Hybrid (MEGAHIT -> Bowtie2 -> Contig Taxonomy)
  # ----------------------------------------
  "assembly_hybrid")
    echo "[INFO] Starting Assembly-Hybrid pipeline..."
    
    # Step A: Assembly
    megahit -1 "$R1" -2 "$R2" -t "$THREADS" -o "$TMP_DIR/megahit" -f
    CONTIGS="$TMP_DIR/megahit/final.contigs.fa"

    if [ ! -s "$CONTIGS" ]; then
        echo "[WARNING] MEGAHIT failed to assemble any contigs. Creating empty output."
        touch "$OUT"
        exit 0
    fi

    # Step B: Map reads to contigs
    bowtie2-build -q "$CONTIGS" "$TMP_DIR/index"
    bowtie2 -x "$TMP_DIR/index" -1 "$R1" -2 "$R2" -p "$THREADS" -S "$TMP_DIR/map.sam" 2>/dev/null

    # Step C: Classify contigs
    kraken2 --db "$DB" --threads "$THREADS" "$CONTIGS" --output "$TMP_DIR/contigs.kraken2" --use-names > /dev/null

    if [ ! -s "$TMP_DIR/contigs.kraken2" ]; then
        echo "[WARNING] Kraken2 did not classify any contigs. Creating empty output."
        touch "$OUT"
        exit 0
    fi

    # Step D: Transfer taxonomy from contigs to reads via Python snippet
    echo "[INFO] Transferring taxonomy from contigs to mapped reads..."
    python3 -c "
import sys
import pandas as pd

sam_file = '$TMP_DIR/map.sam'
kraken_file = '$TMP_DIR/contigs.kraken2'
out_file = '$OUT'

try:
    # Read contig taxonomy (Kraken2 format: status, seq_id, taxid)
    contig_tax = pd.read_csv(kraken_file, sep='\t', header=None, usecols=[1, 2], names=['contig_id', 'taxID'], dtype=str)
    
    # Extract clean taxID
    contig_tax['taxID'] = contig_tax['taxID'].str.extract(r'(?<=taxid )(\d+)|(^\\d+$)')[0].fillna(contig_tax['taxID'].str.extract(r'(?<=taxid )(\d+)|(^\\d+$)')[1])
    tax_dict = dict(zip(contig_tax['contig_id'], contig_tax['taxID']))

    # Parse SAM file and assign taxonomy
    results = []
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'): continue
            parts = line.split('\t')
            read_id = parts[0]
            contig_id = parts[2] # RNAME is 3rd column
            
            # If mapped and contig is classified
            if contig_id != '*' and contig_id in tax_dict and tax_dict[contig_id] != '0':
                results.append(f'{read_id}\t{tax_dict[contig_id]}\n')

    # Write final output
    with open(out_file, 'w') as f:
        f.writelines(results)
except Exception as e:
    print(f'Error merging hybrid taxonomy: {e}', file=sys.stderr)
    sys.exit(1)
"
    ;;

  # ----------------------------------------
  # Fallback
  # ----------------------------------------
  *)
    echo "[ERROR] Unknown tool specified in -p: $TOOL"
    exit 1
    ;;
esac

# Cleanup temporary directory to save disk space
rm -rf "$TMP_DIR"
echo "[SUCCESS] Finished $TOOL wrapper. Results saved to $OUT"

