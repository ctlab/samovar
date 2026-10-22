#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runtime_exports="$("$SCRIPT_DIR/../../bin/samovar" _runtime-env)"
eval "$runtime_exports"

trap 'rm -rf "${TMP_DIR:-}" 2>/dev/null || true' EXIT
# Default threads
THREADS=4
R1="" R2="" DB="" OUT="" TOOL=""

# Consume wrapper options until the native tail. First wrapper values win.
THREADS_SET=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) R1="${R1:-$2}"; shift 2 ;;
    -I) R2="${R2:-$2}"; shift 2 ;;
    -d) DB="${DB:-$2}"; shift 2 ;;
    -o) OUT="${OUT:-$2}"; shift 2 ;;
    -p) TOOL="${TOOL:-$2}"; shift 2 ;;
    -t|--threads)
      if [ "$THREADS_SET" = 0 ]; then THREADS="$2"; THREADS_SET=1; fi
      shift 2 ;;
    --threads=*)
      if [ "$THREADS_SET" = 0 ]; then THREADS="${1#*=}"; THREADS_SET=1; fi
      shift ;;
    --)
      shift
      break
      ;;
    *) break ;;
  esac
done

# Validate required arguments
if [ -z "$R1" ] || [ -z "$OUT" ] || [ -z "$TOOL" ]; then
  echo "[ERROR] Missing required arguments: -i, -o, or -p"
  echo "Usage: ./custom.sh -i R1 [-I R2] -d DB -o OUT -p TOOL_NAME [-t|--threads THREADS]"
  exit 1
fi
# No -I: single-end. Paired tools branch on a non-empty R2 below.

# Remove only adapter-owned options from the native tail. Preserve every
# other argument verbatim, including spaces, order, and literal shell tokens.
EXTRA=()
SEQUENCER=(--sequencer Illumina)
while [ "$#" -gt 0 ]; do
  token="$1"
  key="${token%%=*}"
  owned=0
  boolean=0
  case "$key" in
    -i|-I|-d|-o|-p|-t|--threads) owned=1 ;;
  esac
  case "$TOOL:$key" in
    metax:--db|metax:--dmp-dir|metax:--in-seq|metax:--outprefix) owned=1 ;;
    metax:--is-paired|metax:--no-is-paired) owned=1; boolean=1 ;;
    centrifuge:-x|centrifuge:-1|centrifuge:-2|centrifuge:-U|centrifuge:-S|centrifuge:--report-file) owned=1 ;;
    metax:--sequencer) SEQUENCER=() ;;
  esac
  case "$token" in
    -[iIdopt]?*) owned=1; boolean=1 ;;
  esac
  if [ "$TOOL" = centrifuge ]; then
    case "$token" in -[x12US]?*) owned=1; boolean=1 ;; esac
  fi
  if [ "$owned" = 1 ]; then
    echo "[WARNING] Ignoring extra option owned by SAMOVAR: $token" >&2
    shift
    if [ "$boolean" = 0 ] && [[ "$token" != *=* ]] && [ "$#" -gt 0 ] && [[ "$1" != -* ]]; then
      shift
    fi
  else
    EXTRA+=("$1")
    shift
  fi
done
if [ "$TOOL" = assembly_hybrid ] && [ "${#EXTRA[@]}" -gt 0 ]; then
  echo "[ERROR] assembly_hybrid does not accept a shared native flag tail (assembler and classifier have different CLIs)." >&2
  exit 2
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

# Same lookup as Python: tools.* / tool_envs in the main config, then PATH.
cfg_tool() {
  local name="$1"
  local py="${PYTHON_PATH:-python3}"
  "$py" - "$name" <<'PY'
import os
import sys
from pathlib import Path

from samovar.paths import resolve_executable, command_argv

name = sys.argv[1]
found = command_argv(resolve_executable(name, tool_key=name) or "")
head = found[0] if found else ""
path = Path(head)
if not (path.is_file() and os.access(path, os.X_OK)):
    sys.stderr.write(
        f"[ERROR] {name} is not in the SamovaR config. "
        "Import it with ./install.sh ask or samovar tools import.\n"
    )
    raise SystemExit(127)
print(path)
PY
}

# ==========================================
# TOOL ROUTER
# ==========================================
case "$TOOL" in

  # ----------------------------------------
  # Dummy / constant-taxID classifier (pytest + custom-class smoke test)
  # ----------------------------------------
  "dummy"|"dummy9606"|"constant9606"|"constant"|"random")
    echo "[INFO] Starting constant-taxID dummy classifier (9606)..."
    "$PYTHON_PATH" "$SCRIPT_DIR/constant9606.py" -i "$R1" -I "$R2" -o "$OUT" --taxid 9606 "${EXTRA[@]}"
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
      if [ -s "$R2" ]; then
        ln -sf "$(realpath "$R2")" "$FQ_DIR/R2.fastq.gz"
      fi
    else
      ln -sf "$(realpath "$R1")" "$FQ_DIR/R1.fastq"
      if [ -s "$R2" ]; then
        ln -sf "$(realpath "$R2")" "$FQ_DIR/R2.fastq"
      fi
    fi

    # Run the python script
    "$PYTHON_PATH" "$SCRIPT_DIR/metauto.py" work "$FQ_DIR" "$DB" "${EXTRA[@]}"

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
      if [ -z "$cf1" ]; then
        cf1=$(ls "$DB"/*/*.1.cf 2>/dev/null | head -n 1 || true)
      fi
      if [ -n "$cf1" ]; then
        INDEX="${cf1%.1.cf}"
      fi
    elif [[ "$DB" == *.1.cf ]]; then
      INDEX="${DB%.1.cf}"
    fi
    echo "[INFO] Centrifuge index prefix: ${INDEX}"
    CENTRIFUGE="$(cfg_tool centrifuge)"
    # --mm: the nt index is larger than a job should load into RAM.
    # --report-file: centrifuge otherwise writes centrifuge_report.tsv into CWD.
    # Unpaired CAMI files have no R2: -U, not -1/-2.
    if [ -s "$R2" ]; then
      "$CENTRIFUGE" -x "$INDEX" -1 "$R1" -2 "$R2" -p "$THREADS" --mm \
        -S "$TMP_DIR/centrifuge.out" \
        --report-file "$TMP_DIR/centrifuge_report.tsv" "${EXTRA[@]}"
    else
      "$CENTRIFUGE" -x "$INDEX" -U "$R1" -p "$THREADS" --mm \
        -S "$TMP_DIR/centrifuge.out" \
        --report-file "$TMP_DIR/centrifuge_report.tsv" "${EXTRA[@]}"
    fi
    
    # Parse output: skip header (NR>1), only print classified reads (taxID != 0)
    awk -F'\t' '{if ($3 != 0 && NR>1) print $1 "\t" $3}' "$TMP_DIR/centrifuge.out" > "$OUT"
    ;;

  # ----------------------------------------
  # 4. Metax profiler (classify.txt -> seq, taxID)
  # ----------------------------------------
  "metax")
    echo "[INFO] Starting Metax classification..."
    METAX_BIN="$(cfg_tool metax)"
    # maCMD ships next to the metax binary inside its conda prefix.
    export PATH="$(dirname "$METAX_BIN"):${PATH}"
    DB_JSON="$DB"
    DMP="${METAX_DMP_DIR:-}"
    if [ -d "$DB" ]; then
      if [ -f "$DB/metax_db.json" ]; then
        DB_JSON="$DB/metax_db.json"
      else
        jsons=("$DB"/*.json)
        if [ -f "${jsons[0]}" ]; then
          DB_JSON="${jsons[0]}"
        fi
      fi
      if [ -z "$DMP" ] && [ -f "$DB/nodes.dmp" ]; then
        DMP="$DB"
      elif [ -z "$DMP" ] && [ -f "$DB/taxonomy/nodes.dmp" ]; then
        DMP="$DB/taxonomy"
      fi
    elif [[ "$DB" == *.json ]]; then
      DB_JSON="$DB"
      parent="$(dirname "$DB")"
      if [ -z "$DMP" ] && [ -f "$parent/nodes.dmp" ]; then
        DMP="$parent"
      elif [ -z "$DMP" ] && [ -f "$parent/taxonomy/nodes.dmp" ]; then
        DMP="$parent/taxonomy"
      fi
    fi
    if [ -z "$DMP" ] || [ ! -f "$DMP/nodes.dmp" ]; then
      py="${PYTHON_PATH:-python3}"
      resolved="$("$py" - <<'PY'
from samovar.taxdump import metax_dmp_dir
print(metax_dmp_dir() or "")
PY
)"
      if [ -n "$resolved" ] && [ -f "$resolved/nodes.dmp" ]; then
        DMP="$resolved"
      fi
    fi
    if [ -z "$DMP" ] || [ ! -f "$DMP/nodes.dmp" ]; then
      echo "[ERROR] Metax needs a taxonomy dump (--dmp-dir). Set genomes.taxdump / databases.taxdump in the SamovaR config, METAX_DMP_DIR, or place nodes.dmp next to the database." >&2
      exit 2
    fi
    # Point the Metax DB at the configured dump when it does not have its own.
    parent=""
    if [ -d "$DB" ]; then
      parent="$DB"
    elif [[ "$DB" == *.json ]]; then
      parent="$(dirname "$DB")"
    fi
    if [ -n "$parent" ] && [ ! -e "$parent/taxonomy" ] && [ ! -e "$parent/nodes.dmp" ] && [ -w "$parent" ]; then
      ln -s "$DMP" "$parent/taxonomy" || true
    fi
    echo "[INFO] Metax db=${DB_JSON} dmp=${DMP}"
    PREFIX="$TMP_DIR/metax"
    if [ -s "$R2" ]; then
      "$METAX_BIN" profile --db "$DB_JSON" --dmp-dir "$DMP" \
        --in-seq "${R1},${R2}" --is-paired "${SEQUENCER[@]}" \
        --outprefix "$PREFIX" --threads "$THREADS" "${EXTRA[@]}"
    else
      "$METAX_BIN" profile --db "$DB_JSON" --dmp-dir "$DMP" \
        --in-seq "$R1" "${SEQUENCER[@]}" \
        --outprefix "$PREFIX" --threads "$THREADS" "${EXTRA[@]}"
    fi
    CLASSIFY="${PREFIX}.classify.txt"
    if [ ! -s "$CLASSIFY" ]; then
      echo "[WARNING] Metax wrote no classify.txt (fractional indexes omit it)."
      touch "$OUT"
    else
      awk -F'\t' 'NR==1 && ($2 !~ /^[0-9]+$/ || $1 ~ /^(read|seq|query)$/) {next} $2 ~ /^[0-9]+$/ && $2 != "0" {print $1 "\t" $2}' \
        "$CLASSIFY" > "$OUT"
    fi
    ;;

  # ----------------------------------------
  # 3. Assembly + Mapping Hybrid (MEGAHIT -> Bowtie2 -> Contig Taxonomy)
  # ----------------------------------------
  "assembly_hybrid")
    echo "[INFO] Starting Assembly-Hybrid pipeline..."
    
    KRAKEN2="$(cfg_tool kraken2)"
    CONTIGS="$TMP_DIR/contigs.fasta"
    "${PYTHON_PATH:-python3}" - "$R1" "${R2:-}" "$CONTIGS" "$THREADS" "$TMP_DIR/map.bam" <<'PYCODE'
import sys
from pathlib import Path
from samovar.assembly_profiling import run_assembler, run_minimap2, which_tool, _run
r1, r2, contigs, threads, bam = sys.argv[1:]
run_assembler(r1, r2, contigs, threads=int(threads), name="spades")
if Path(contigs).stat().st_size:
    run_minimap2(r1, r2, contigs, bam, threads=int(threads))
    _run([which_tool("samtools"), "view", "-h", "-o", str(Path(bam).with_suffix(".sam")), bam])
PYCODE
    if [ ! -s "$CONTIGS" ]; then
        echo "[WARNING] SPAdes produced no contigs. Creating empty output."
        touch "$OUT"
        exit 0
    fi

    # Step C: Classify contigs
    "$KRAKEN2" --db "$DB" --threads "$THREADS" "$CONTIGS" --output "$TMP_DIR/contigs.kraken2" --use-names > /dev/null

    if [ ! -s "$TMP_DIR/contigs.kraken2" ]; then
        echo "[WARNING] Kraken2 did not classify any contigs. Creating empty output."
        touch "$OUT"
        exit 0
    fi

    # Step D: Transfer taxonomy from contigs to reads via Python snippet
    echo "[INFO] Transferring taxonomy from contigs to mapped reads..."
    "${PYTHON_PATH:-python3}" -c "
import sys
import pandas as pd

sam_file, kraken_file, out_file = sys.argv[1:]

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
" "$TMP_DIR/map.sam" "$TMP_DIR/contigs.kraken2" "$OUT"
    ;;

  # ----------------------------------------
  # Fallback
  # ----------------------------------------
  *)
    echo "[ERROR] Unknown tool specified in -p: $TOOL"
    exit 1
    ;;
esac

# Cleanup temporary directory to save disk space (NFS .nfs* files may be busy)
rm -rf "$TMP_DIR" || true
echo "[SUCCESS] Finished $TOOL wrapper. Results saved to $OUT"

