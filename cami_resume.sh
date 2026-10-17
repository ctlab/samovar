#!/usr/bin/env bash
# Stage finished CAMI initial profiles, then regenerate one 100k-read sample
# and reprofile it. Visualizations and the full initial table are skipped.
#
#   THREADS=200 bash cami_resume.sh
set -euo pipefail

CAMI="${CAMI:-/mnt/san03/gentech/a.taldaev/cami}"
OUT="${OUT:-$CAMI/samovar/resume}"
READS="${READS:-$CAMI/raw_sr}"
SR="${SR:-$CAMI/sr}"
N_KEEP="${N_KEEP:-500}"
TOP_TAXA="${TOP_TAXA:-150}"
MAX_GENOME_MB="${MAX_GENOME_MB:-100}"
N_READS="${N_READS:-100000}"
THREADS="${THREADS:-64}"

METAX_DB="${METAX_DB:-/nvme/a.taldaev/metax-abfph_cc-v_ccs-refseq_202509}"
KAIJU_DB="${KAIJU_DB:-/mnt/san02/gentech/a.taldaev/db/kaijudb/refseq}"
KAIJU_FMI="${KAIJU_FMI:-$KAIJU_DB/kaiju_db_refseq.fmi}"
KAIJU_NODES="${KAIJU_NODES:-/mnt/san02/gentech/a.taldaev/db/kaijudb/nodes.dmp}"
KRAKEN2_DB="${KRAKEN2_DB:-/dev/shm/k2_core_nt}"
CENTRIFUGE_DB="${CENTRIFUGE_DB:-/nvme/a.taldaev/centrifuge_nt}"
GTDBTK_DB="${GTDBTK_DB:-/nvme/a.taldaev/GTDB_r232}"

# This cluster has no GPU. A previous prepare baked --gpu into samovar.sh and
# gpu: true into the reprofile config; both would demand cuML.
drop_pipeline_gpu() {
    local dir="$OUT/.log"
    local script yaml
    [ -d "$dir" ] || return 0
    shopt -s nullglob
    for script in "$dir"/samovar.sh "$dir"/samovar_v*.sh; do
        [ -f "$script" ] || continue
        sed -i -E '/^[[:space:]]*--gpu[[:space:]]*\\?[[:space:]]*$/d' "$script"
    done
    shopt -u nullglob
    for yaml in "$dir"/configs/*.yaml "$OUT"/.hydra/*.yaml; do
        [ -f "$yaml" ] || continue
        sed -i -E "s/^([[:space:]]*gpu:[[:space:]]*)['\"]?[Tt]rue['\"]?[[:space:]]*$/\\1false/" "$yaml"
    done
    echo "[resume] RandomForest on CPU"
}

# Annotation of the regenerated reads finished. Continue at the table join,
# or at reprofile once that join is done.
CFG="$OUT/.log/configs/config_reannotate.yaml"
ASSEMBLY_OUT="$OUT/regenerated_reports/pool_concat_anvio.assembly.out"
if [ -f "$OUT/.log/checkpoints/annotate_regenerated.done" ] && [ -s "$ASSEMBLY_OUT" ]; then
    PY="${SAMOVAR_PYTHON:-}"
    if [ -z "$PY" ] && [ -x /nvme/a.taldaev/miniconda3/envs/samovar/bin/python3 ]; then
        PY=/nvme/a.taldaev/miniconda3/envs/samovar/bin/python3
    fi
    if [ -z "$PY" ]; then
        PY="$(command -v python3 || true)"
    fi
    export OUT PY
    "$PY" - <<'PY'
import os
import shutil
from pathlib import Path

out = Path(os.environ["OUT"])
reports = out / "initial_reports"
initial = out / "initial_annotations"
regen = out / "regenerated_annotations"
initial.mkdir(parents=True, exist_ok=True)
regen.mkdir(parents=True, exist_ok=True)

def csvs(folder):
    return [p for p in folder.glob("*.csv") if p.is_file() and not p.name.startswith(".")]

if not csvs(initial):
    outs = list(reports.glob("*.out")) if reports.is_dir() else []
    if not outs:
        raise SystemExit(f"No initial reports in {reports}")
    from samovar.combine_tables import combine_with_cpp
    print(f"[resume] combining {len(outs)} initial reports into {initial}")
    combine_with_cpp(str(reports), str(initial), 1, 500000)

combined = regen / "combined_annotation_table.csv"
parts = sorted(p for p in regen.glob("*.annotation.csv") if p.is_file())
if not combined.is_file() or combined.stat().st_size == 0:
    if len(parts) == 1:
        shutil.copyfile(parts[0], combined)
        print(f"[resume] {combined.name} <- {parts[0].name}")
    elif len(parts) > 1:
        import pandas as pd
        pd.concat([pd.read_csv(p) for p in parts], ignore_index=True).to_csv(combined, index=False)
        print(f"[resume] concatenated {len(parts)} regenerated tables -> {combined.name}")
    else:
        raise SystemExit(f"No regenerated annotation table in {regen}")

# Same tax/feature columns on every initial sample as on the regenerated table.
import pandas as pd
header = list(pd.read_csv(combined, nrows=0).columns)
contract = [c for c in header if str(c).startswith("taxID_") or str(c).startswith("feat_")]
for path in csvs(initial):
    if path.name.startswith("combined_annotation_table"):
        continue
    frame = pd.read_csv(path)
    missing = [c for c in contract if c not in frame.columns]
    if not missing:
        continue
    for col in missing:
        frame[col] = 0 if str(col).startswith("taxID_") else ""
    frame.to_csv(path, index=False)
    print(f"[resume] {path.name}: added {', '.join(missing)}")
PY
    # Plots and MultiQC were stamped done up front. Run them after reprofile.
    rm -f "$OUT/.log/checkpoints/viz_regenerated.done" \
          "$OUT/.log/checkpoints/viz_reprofiled.done"
    export SAMOVAR_MULTIQC=1
    if [ -f "$OUT/.log/checkpoints/combine_regenerated.done" ]; then
        echo "[resume] combine_regenerated is done; starting at reprofile"
    else
        echo "[resume] annotate_regenerated is done; starting at combine_regenerated"
    fi
    drop_pipeline_gpu
    set +e
    samovar exec --directory "$OUT"
    status=$?
    set -e
    exit "$status"
fi

# Regeneration already finished: only correct the databases and repeat annotation.
if [ -f "$CFG" ] && compgen -G "$OUT/regenerated_trimmed/*_R1.fastq" >/dev/null; then
    for piece in hash.k2d opts.k2d taxo.k2d; do
        if [ ! -e "$KRAKEN2_DB/$piece" ]; then
            echo "Kraken2 database $KRAKEN2_DB has no $piece." >&2
            echo "kraken2 then exits 2 before reading any FASTQ. Restore that directory in /dev/shm and rerun." >&2
            exit 1
        fi
    done
    if ! compgen -G "$CENTRIFUGE_DB/*.1.cf" >/dev/null && ! compgen -G "$CENTRIFUGE_DB/*/*.1.cf" >/dev/null; then
        echo "Centrifuge index $CENTRIFUGE_DB has no *.1.cf" >&2
        exit 1
    fi
    for need in "$KAIJU_FMI" "$KAIJU_NODES" "$METAX_DB" "$GTDBTK_DB"; do
        if [ ! -e "$need" ]; then
            echo "Missing database path: $need" >&2
            exit 1
        fi
    done

    # The installed annotator calls libexec/kraken2 --use-names. The manual run
    # that produced ./k2/*.out.gz used bin/kraken2 without --use-names. A /dev/shm
    # copy often has only the *.k2d files, and --use-names then exits 2 at once.
    K2_WRAP="$OUT/.log/kraken2.sh"
    mkdir -p "$OUT/.log"
    cat > "$K2_WRAP" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
real="${KRAKEN2_BIN:-}"
if [ -z "$real" ]; then
    real="$(command -v kraken2 || true)"
fi
if [ -z "$real" ]; then
    echo "kraken2 is not on PATH" >&2
    exit 127
fi
# Prefer the conda bin wrapper over libexec/kraken2.
if [[ "$real" == */libexec/kraken2 ]]; then
    env_root="${real%/share/*}"
    if [ -x "$env_root/bin/kraken2" ]; then
        real="$env_root/bin/kraken2"
    fi
fi
db=""
prev=""
for arg in "$@"; do
    if [ "$prev" = "--db" ]; then
        db="$arg"
    fi
    prev="$arg"
done
drop_names=0
if [ -n "$db" ] && [ ! -f "$db/taxonomy/names.dmp" ]; then
    drop_names=1
fi
args=()
out=""
prev=""
for arg in "$@"; do
    if [ "$drop_names" = 1 ] && [ "$arg" = "--use-names" ]; then
        prev="$arg"
        continue
    fi
    if [ "$prev" = "--output" ]; then
        out="$arg"
    fi
    args+=("$arg")
    prev="$arg"
done
"$real" "${args[@]}"
if [ "$drop_names" = 1 ] && [ -n "$out" ] && [ -f "$out" ]; then
    # Installed parser only keeps a taxid written as "taxid N".
    awk -F '\t' 'BEGIN{OFS="\t"} $3 ~ /^[0-9]+$/ {$3=$3" (taxid "$3")"} {print}' \
        "$out" > "$out.named" && mv "$out.named" "$out"
fi
EOF
    chmod +x "$K2_WRAP"
    if [ -x /nvme/a.taldaev/miniconda3/envs/k2/bin/kraken2 ]; then
        export KRAKEN2_BIN=/nvme/a.taldaev/miniconda3/envs/k2/bin/kraken2
    fi

    # map-initial still ran megahit, then loaded every MAG FASTA into one Python
    # string before minimap2. That is the assembly job that exits 1. Map the
    # regenerated reads onto the staged MAG fastas directly.
    # Executables come from the SamovaR main config (tools.minimap2 / tools.samtools),
    # the same lookup the assembly annotator uses.
    PY="${SAMOVAR_PYTHON:-}"
    if [ -z "$PY" ] && [ -x /nvme/a.taldaev/miniconda3/envs/samovar/bin/python3 ]; then
        PY=/nvme/a.taldaev/miniconda3/envs/samovar/bin/python3
    fi
    if [ -z "$PY" ]; then
        PY="$(command -v python3 || true)"
    fi
    if [ -z "$PY" ]; then
        echo "python3 from the samovar env is required to read the main config." >&2
        exit 1
    fi
    tool_bins="$("$PY" - <<'PY'
from samovar.assembly_profiling import which_tool
for name in ("minimap2", "samtools"):
    print(which_tool(name))
PY
)" || {
        echo "minimap2 or samtools is missing from the SamovaR main config." >&2
        exit 1
    }
    MINIMAP2="$(printf '%s\n' "$tool_bins" | sed -n '1p')"
    SAMTOOLS="$(printf '%s\n' "$tool_bins" | sed -n '2p')"
    if [ ! -x "$MINIMAP2" ] || [ ! -x "$SAMTOOLS" ]; then
        echo "Main config paths are not executable: minimap2=$MINIMAP2 samtools=$SAMTOOLS" >&2
        exit 1
    fi
    MAP_WRAP="$OUT/.log/map_mags.sh"
    cat > "$MAP_WRAP" <<EOF
#!/usr/bin/env bash
set -euo pipefail
minimap2="$MINIMAP2"
samtools="$SAMTOOLS"
mags="\${SAMOVAR_MAG_FASTAS:-$OUT/MAG/fastas}"
r1="" out="" threads=4 work=""
args=("\$@")
i=0
while [ "\$i" -lt "\${#args[@]}" ]; do
    a="\${args[\$i]}"
    nxt="\${args[\$((i+1))]:-}"
    case "\$a" in
        -i) r1="\$nxt"; i=\$((i+2)) ;;
        -o) out="\$nxt"; i=\$((i+2)) ;;
        -t) threads="\$nxt"; i=\$((i+2)) ;;
        --work-dir) work="\$nxt"; i=\$((i+2)) ;;
        *) i=\$((i+1)) ;;
    esac
done
if [ -z "\$r1" ] || [ -z "\$out" ]; then
    echo "map_mags.sh needs -i reads and -o table" >&2
    exit 2
fi
if [ -z "\$work" ]; then
    work="\${out%.out}"
fi
mkdir -p "\$work" "\$(dirname "\$out")"
shopt -s nullglob
fastas=("\$mags"/*.fa)
if [ "\${#fastas[@]}" -eq 0 ]; then
    echo "No MAG FASTA in \$mags" >&2
    exit 1
fi
ref="\$work/mags.fa"
: > "\$ref"
for fa in "\${fastas[@]}"; do
    mag="\$(basename "\$fa" .fa)"
    awk -v mag="\$mag" '
        /^>/ {
            contig = substr(\$1, 2)
            if (index(contig, mag "~") == 1) print ">" contig
            else print ">" mag "~" contig
            next
        }
        { print }
    ' "\$fa" >> "\$ref"
done
bam="\$work/reads_to_mag.bam"
# A reference this large is indexed in parts. Without --split-prefix
# minimap2 omits @SQ lines and samtools sort aborts on the truncated SAM.
"\$minimap2" -ax sr -t "\$threads" --split-prefix "\$work/mm2split" "\$ref" "\$r1" \\
    | "\$samtools" sort -@ "\$threads" -o "\$bam"
"\$samtools" index "\$bam" || true
{
    printf 'seq\\ttaxID\\tMAP_MAG_ID\\n'
    "\$samtools" view -F 2308 "\$bam" \\
        | awk -F '\\t' 'BEGIN{OFS="\\t"} NF>=3 && \$3!="*" { mag=\$3; sub(/~.*/, "", mag); print \$1, "0", mag }'
} > "\$out"
EOF
    chmod +x "$MAP_WRAP"
    export SAMOVAR_MAG_FASTAS="$OUT/MAG/fastas"

    export OUT CFG METAX_DB KAIJU_FMI KAIJU_NODES KRAKEN2_DB CENTRIFUGE_DB GTDBTK_DB K2_WRAP MAP_WRAP
    python3 - <<'PY'
import os
from pathlib import Path
import yaml
cfg_path = Path(os.environ["CFG"])
data = yaml.safe_load(cfg_path.read_text()) or {}
wanted = {
    "metax": os.environ["METAX_DB"],
    "kaiju": os.environ["KAIJU_FMI"],
    "kraken2": os.environ["KRAKEN2_DB"],
    "centrifuge": os.environ["CENTRIFUGE_DB"],
    "anvio": os.environ["GTDBTK_DB"],
}
custom = ""
for row in data.get("run_config") or []:
    name = str(row.get("run_name") or row.get("type") or "")
    if name not in wanted:
        continue
    row["db_path"] = wanted[name]
    if name == "kaiju":
        row["db_nodes"] = os.environ["KAIJU_NODES"]
    if name == "kraken2":
        row["cmd"] = os.environ["K2_WRAP"]
        extra = str(row.get("extra") or "")
        if "--memory-mapping" not in extra:
            row["extra"] = (extra + " --memory-mapping").strip()
    if name == "centrifuge":
        for token in str(row.get("cmd") or "").split():
            if token.endswith("custom.sh"):
                custom = token
    if name == "anvio":
        row["cmd"] = os.environ["MAP_WRAP"]
    print(f"{name}: {row['db_path']}")
cfg_path.write_text(yaml.safe_dump(data, sort_keys=False))
if custom:
    Path(os.environ["OUT"] + "/.log/centrifuge_custom.sh").write_text(custom)
PY
    CUSTOM=""
    if [ -f "$OUT/.log/centrifuge_custom.sh" ]; then
        CUSTOM="$(cat "$OUT/.log/centrifuge_custom.sh")"
    fi
    if [ -n "$CUSTOM" ] && [ -f "$CUSTOM" ] && ! grep -q -- '--mm' "$CUSTOM"; then
        python3 - "$CUSTOM" <<'PY'
import sys
from pathlib import Path
path = Path(sys.argv[1])
text = path.read_text()
old = '"$CENTRIFUGE"'
if "--mm" in text or old not in text:
    raise SystemExit(0)
path.write_text(text.replace('-p "$THREADS" \\', '-p "$THREADS" --mm \\'))
print(f"added --mm to {path}")
PY
    fi
    assembly_out="$OUT/regenerated_reports/pool_concat_anvio.assembly.out"
    if [ ! -s "$assembly_out" ]; then
        rm -f "$assembly_out"
        rm -rf "$OUT/regenerated_reports/pool_concat_anvio.assembly"
        rm -f "$OUT/.log/checkpoints/annotate_regenerated.done"
    fi
    export SAMOVAR_MAG_MODE=map-initial
    drop_pipeline_gpu
    set +e
    samovar exec --directory "$OUT"
    status=$?
    set -e
    if [ "$status" != 0 ]; then
        alog="$OUT/regenerated_reports/logs/pool_concat_anvio.assembly.log"
        if [ ! -s "$OUT/regenerated_reports/pool_concat_anvio.assembly.out" ] && [ -s "$alog" ]; then
            echo "---- assembly ----" >&2
            cat "$alog" >&2
        fi
        exit "$status"
    fi
    exit 0
fi

if [ ! -f "$SR" ]; then
    echo "Missing sample list: $SR" >&2
    exit 1
fi

mkdir -p \
    "$OUT/reads" "$OUT/initial" "$OUT/initial_trimmed" \
    "$OUT/initial_reports" "$OUT/initial_abundance" \
    "$OUT/external" "$OUT/MAG/fastas" \
    "$OUT/.log/checkpoints"

hardlink() {
    local src="$1" dest="$2"
    mkdir -p "$(dirname "$dest")"
    if [ -e "$dest" ] || [ -L "$dest" ]; then
        return 0
    fi
    cp -l "$src" "$dest" 2>/dev/null || ln -s "$src" "$dest"
}

while read -r id; do
    [ -n "$id" ] || continue
    src="$READS/${id}.fq.gz"
    if [ ! -f "$src" ]; then
        echo "Missing reads: $src" >&2
        exit 1
    fi
    for dest in reads initial initial_trimmed; do
        target="$OUT/$dest/${id}_R1.fq.gz"
        if [ ! -e "$target" ] && [ ! -L "$target" ]; then
            ln -s "$src" "$target"
        fi
    done
    [ -f "$CAMI/kaiju/${id}.out" ] && hardlink "$CAMI/kaiju/${id}.out" "$OUT/external/kaiju/${id}.out"
    if [ -f "$CAMI/k2/${id}.out.gz" ]; then
        hardlink "$CAMI/k2/${id}.out.gz" "$OUT/external/k2/${id}.out.gz"
        hardlink "$CAMI/k2/${id}.report" "$OUT/external/k2/${id}.report"
    fi
    if [ -f "$CAMI/centrifuge/${id}.out.gz" ]; then
        hardlink "$CAMI/centrifuge/${id}.out.gz" "$OUT/external/centrifuge/${id}.out.gz"
        hardlink "$CAMI/centrifuge/${id}.report" "$OUT/external/centrifuge/${id}.report"
    fi
    [ -f "$CAMI/metax/${id}.classify.txt" ] && hardlink "$CAMI/metax/${id}.classify.txt" "$OUT/external/metax/${id}.classify.txt"
    [ -f "$CAMI/metax/metax_host/${id}.classify.txt" ] && hardlink "$CAMI/metax/metax_host/${id}.classify.txt" "$OUT/external/metax_host/${id}.classify.txt"
    mag=""
    if [ -f "$CAMI/anvi-fasta/${id}.contigs.fa" ]; then
        mag="$CAMI/anvi-fasta/${id}.contigs.fa"
    elif [ -f "$CAMI/megahit/${id}/${id}.contigs.fa" ]; then
        mag="$CAMI/megahit/${id}/${id}.contigs.fa"
    fi
    [ -n "$mag" ] && hardlink "$mag" "$OUT/MAG/fastas/${id}.fa"
    [ -f "$CAMI/mapping/${id}.bam" ] && hardlink "$CAMI/mapping/${id}.bam" "$OUT/external/mapping/${id}.bam"
done < "$SR"

if [ ! -f "$OUT/MAG/catalog.tsv" ]; then
    printf 'mag_id\tn\tcounts\tsource\n' > "$OUT/MAG/catalog.tsv"
fi

export CAMI OUT SR N_KEEP TOP_TAXA
python3 - <<'PY'
import gzip
import os
from collections import Counter
from pathlib import Path

cami = Path(os.environ["CAMI"])
out = Path(os.environ["OUT"])
sr = Path(os.environ["SR"])
n_keep = int(os.environ["N_KEEP"])
top = int(os.environ["TOP_TAXA"])
ids = [line.strip() for line in sr.read_text().splitlines() if line.strip()]
reports = out / "initial_reports"
abundance = out / "initial_abundance"
reports.mkdir(parents=True, exist_ok=True)
abundance.mkdir(parents=True, exist_ok=True)


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def take(pairs, limit, scan_cap=200000):
    kept = []
    seen = 0
    for seq, tax in pairs:
        seen += 1
        if seen > scan_cap:
            break
        if not tax or tax == "0":
            continue
        kept.append((seq, tax))
        if len(kept) >= limit:
            break
    return kept


def kaiju_pairs(path: Path):
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != "C":
                continue
            yield parts[1], parts[2]


def kraken_pairs(path: Path):
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != "C":
                continue
            taxa = parts[2].strip()
            if taxa.isdigit():
                yield parts[1], taxa
                continue
            mark = taxa.rfind("taxid ")
            if mark < 0:
                continue
            tax = taxa[mark + 6 :].split(")", 1)[0].strip()
            if tax.isdigit():
                yield parts[1], tax


def centrifuge_pairs(path: Path):
    with open_text(path) as handle:
        first = True
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if first:
                first = False
                if parts and parts[0].lower() in {"readid", "read"}:
                    continue
            if len(parts) < 3:
                continue
            yield parts[0], parts[2]


def metax_pairs(path: Path):
    with open_text(path) as handle:
        first = True
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if first:
                first = False
                if not parts[1].isdigit():
                    continue
            yield parts[0], parts[1]


def first_existing(candidates):
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            return path
    return None


sources = {
    "kaiju": lambda sample: first_existing([cami / "kaiju" / f"{sample}.out"]),
    "kraken2": lambda sample: first_existing(
        [cami / "k2" / f"{sample}.out.gz", cami / "k2" / f"{sample}.out"]
    ),
    "centrifuge": lambda sample: first_existing(
        [cami / "centrifuge" / f"{sample}.out.gz", cami / "centrifuge" / f"{sample}.out"]
    ),
    "metax": lambda sample: first_existing(
        [
            cami / "metax" / "metax_host" / f"{sample}.classify.txt",
            cami / "metax" / f"{sample}.classify.txt",
        ]
    ),
}
readers = {
    "kaiju": kaiju_pairs,
    "kraken2": kraken_pairs,
    "centrifuge": centrifuge_pairs,
    "metax": metax_pairs,
}
names = {
    "kaiju": "{sample}_kaiju.kaiju.out",
    "kraken2": "{sample}_kraken2.kraken2.out",
    "centrifuge": "{sample}_centrifuge.custom_centrifuge.out",
    "metax": "{sample}_metax.custom_metax.out",
}

for tool, locate in sources.items():
    counts = Counter()
    used = 0
    for sample in ids:
        src = locate(sample)
        if src is None:
            print(f"[skip] {tool} {sample}: no file")
            continue
        rows = take(readers[tool](src), n_keep)
        dest = reports / names[tool].format(sample=sample)
        with dest.open("w", encoding="utf-8") as handle:
            handle.write("seq\ttaxID\n")
            for seq, tax in rows:
                handle.write(f"{seq}\t{tax}\n")
                counts[tax] += 1
        used += 1
    table = abundance / f"{tool}.csv"
    with table.open("w", encoding="utf-8") as handle:
        handle.write("taxid,N_pool\n")
        for tax, n in counts.most_common(top):
            handle.write(f"{tax},{n}\n")
    print(f"[abundance] {tool}: {used} samples, {min(top, len(counts))} taxids -> {table}")
PY

stamp() {
    printf '%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/.log/checkpoints/$1.done"
}
for step in \
    setup_reads qc_initial annotate_initial combine_initial viz_initial \
    abundance_tables viz_regenerated viz_reprofiled
do
    stamp "$step"
done

samovar prepare \
    --input_dir "$OUT/reads" \
    --directory "$OUT" \
    --cores "$THREADS" \
    --startpoint regenerate_tables \
    --endpoint reprofile \
    --max-genomes "$TOP_TAXA" \
    --max-genome-mb "$MAX_GENOME_MB" \
    --table_reads_generator concat \
    --N 1 \
    --N_reads "$N_READS" \
    --reads_generator iss \
    --metax "metax ${METAX_DB}" \
    --kaiju "kaiju ${KAIJU_DB}" \
    --kraken2 "kraken2 ${KRAKEN2_DB} --memory-mapping" \
    --centrifuge "centrifuge ${CENTRIFUGE_DB}" \
    --anvio "assembly_profiling ${GTDBTK_DB}" \
    --assembler megahit \
    --binner anvio \
    --mag-taxonomy identity \
    --mag-mode map-initial \
    --to-cami \
    --taxonomy gtdb

export SAMOVAR_MAG_MODE=map-initial
drop_pipeline_gpu
samovar exec --directory "$OUT"
