# Запуск SAMOVAR на CAMI (другой кластер)

Три независимых скрипта. В каждом пути этого кластера задаются переменными в начале. Данные CAMI III Toy Human Gut лежат двумя наборами с одним и тем же рисунком имён:

- Illumina, `short/`: `sample_{0..19}_1.fastq` и `sample_{0..19}_2.fastq` (пары 2×150).
- ONT, `long/`: те же `sample_{N}_1.fastq` / `_2.fastq`.

SAMOVAR берёт имя образца как всё, что стоит перед `_1` / `_R1`. Без постфикса оба запуска называются `sample_0`, и `merge` их не склеит. Скрипты складывают в рабочую папку ссылки `sample_{N}_illumina_1.fastq` и `sample_{N}_ont_1.fastq`. Если на кластере постфикс уже есть в имени, оставьте `STAGE_READS=0` и укажите `READS` на готовую папку.

Общее для обоих профилирований:

- один метод регенерации abundance table, `concat`: профили всех аннотаторов сливаются в одну таблицу, а не в таблицу на каждый профилировщик;
- `--N 1` — один синтетический образец в конфиге регенерации;
- `--N_reads 10000000`;
- не больше 100 геномов и 500 Мб на новый скачанный геном;
- модель `ensemble`, её же третий скрипт обучает заново на объединённой таблице.

`concat` оставляет в этой одной таблице колонку на каждый входной образец. `--N 1` ограничивает число синтетических профилей у методов вроде bootstrap; у `concat` колонки входных образцов сохраняются.

## 1. Только Illumina

MetaX, Kaiju, Kraken 2, Centrifuge, anvi'o (MEGAHIT, картирование на стартовые MAG, фича `MAP_MAG_ID`, GTDB-Tk), k-mer профиль k=4, GC.

```bash
#!/usr/bin/env bash
set -euo pipefail

# --- пути этого кластера ---
SAMOVAR=/path/to/samovar
READS_SRC=/path/to/cami/short          # sample_N_1.fastq + sample_N_2.fastq
STAGE_READS=1                          # 0, если в именах уже есть _illumina_
READS=/path/to/work/illumina_reads
OUT=/path/to/work/samovar_illumina

METAX_DB=/path/to/metax_db
KAIJU_DB=/path/to/kaiju_db
KRAKEN2_DB=/path/to/kraken2_db
CENTRIFUGE_DB=/path/to/centrifuge_db
GTDBTK_DB=/path/to/gtdbtk_db

CORES=32
MAX_GENOMES=100
MAX_GENOME_MB=500
N_READS=10000000
NCBI_EMAIL=you@example.com
# Раскомментируйте, чтобы не качать хозяина (в toy gut это taxid 9606).
# GENOME_SKIP_LIST=9606
GENOME_SKIP_LIST=""

export NCBI_EMAIL
export PATH="${SAMOVAR}/bin:${PATH}"
export PYTHONPATH="${SAMOVAR}/src${PYTHONPATH:+:${PYTHONPATH}}"
export GTDBTK_DATA_PATH="${GTDBTK_DB}"

if [[ "${STAGE_READS}" == "1" ]]; then
  mkdir -p "${READS}"
  shopt -s nullglob
  for src in "${READS_SRC}"/*; do
    base="$(basename "${src}")"
    [[ -f "${src}" ]] || continue
    case "${base}" in
      *_1.fastq|*_1.fastq.gz|*_1.fq|*_1.fq.gz|*_R1.fastq|*_R1.fastq.gz|*_R1.fq|*_R1.fq.gz)
        mate="${base##*_}"
        stem="${base%_${mate}}"
        if [[ "${stem}" == *_illumina ]]; then
          ln -sfn "${src}" "${READS}/${base}"
        else
          ln -sfn "${src}" "${READS}/${stem}_illumina_${mate}"
        fi
        ;;
      *_2.fastq|*_2.fastq.gz|*_2.fq|*_2.fq.gz|*_R2.fastq|*_R2.fastq.gz|*_R2.fq|*_R2.fq.gz)
        mate="${base##*_}"
        stem="${base%_${mate}}"
        if [[ "${stem}" == *_illumina ]]; then
          ln -sfn "${src}" "${READS}/${base}"
        else
          ln -sfn "${src}" "${READS}/${stem}_illumina_${mate}"
        fi
        ;;
    esac
  done
fi

skip_args=()
if [[ -n "${GENOME_SKIP_LIST}" ]]; then
  skip_args=(--genome-skip-list "${GENOME_SKIP_LIST}")
fi

samovar prepare \
  --input_dir "${READS}" \
  --output_dir "${OUT}" \
  --cores "${CORES}" \
  --max-genomes "${MAX_GENOMES}" \
  --max-genome-mb "${MAX_GENOME_MB}" \
  "${skip_args[@]}" \
  --N 1 \
  --N_reads "${N_READS}" \
  --table_reads_generator concat \
  --reads_generator iss \
  --qc fastp \
  --reprofiler ensemble \
  --metax "metax ${METAX_DB}" \
  --kaiju "kaiju ${KAIJU_DB}" \
  --kraken2 "kraken2 ${KRAKEN2_DB}" \
  --centrifuge "centrifuge ${CENTRIFUGE_DB}" \
  --anvio "assembly_profiling ${GTDBTK_DB}" \
  --assembler megahit \
  --binner anvio \
  --mag-taxonomy gtdbtk \
  --mag-mode map-initial \
  --kmer "kmer_encoder . -k 4" \
  --gc "gc ."

samovar exec --output_dir "${OUT}"
```

## 2. Только ONT

MetaX, Kaiju, Kraken 2, длина рида, k-mer профиль k=4. Регенерация та же: один `concat`, `--N 1`, 10 млн прочтений. Симулятор регенерации — NanoSim, QC — Chopper.

```bash
#!/usr/bin/env bash
set -euo pipefail

# --- пути этого кластера ---
SAMOVAR=/path/to/samovar
READS_SRC=/path/to/cami/long
STAGE_READS=1
READS=/path/to/work/ont_reads
OUT=/path/to/work/samovar_ont

METAX_DB=/path/to/metax_db
KAIJU_DB=/path/to/kaiju_db
KRAKEN2_DB=/path/to/kraken2_db
# Каталог модели NanoSim; пусто — встроенный вызов без --model.
NANOSIM_MODEL=""

CORES=32
MAX_GENOMES=100
MAX_GENOME_MB=500
N_READS=10000000
NCBI_EMAIL=you@example.com
GENOME_SKIP_LIST=""

export NCBI_EMAIL
export PATH="${SAMOVAR}/bin:${PATH}"
export PYTHONPATH="${SAMOVAR}/src${PYTHONPATH:+:${PYTHONPATH}}"

if [[ "${STAGE_READS}" == "1" ]]; then
  mkdir -p "${READS}"
  shopt -s nullglob
  for src in "${READS_SRC}"/*; do
    base="$(basename "${src}")"
    [[ -f "${src}" ]] || continue
    case "${base}" in
      *_1.fastq|*_1.fastq.gz|*_1.fq|*_1.fq.gz|*_R1.fastq|*_R1.fastq.gz|*_R1.fq|*_R1.fq.gz|\
      *_2.fastq|*_2.fastq.gz|*_2.fq|*_2.fq.gz|*_R2.fastq|*_R2.fastq.gz|*_R2.fq|*_R2.fq.gz)
        mate="${base##*_}"
        stem="${base%_${mate}}"
        if [[ "${stem}" == *_ont ]]; then
          ln -sfn "${src}" "${READS}/${base}"
        else
          ln -sfn "${src}" "${READS}/${stem}_ont_${mate}"
        fi
        ;;
    esac
  done
fi

skip_args=()
if [[ -n "${GENOME_SKIP_LIST}" ]]; then
  skip_args=(--genome-skip-list "${GENOME_SKIP_LIST}")
fi
nanosim_args=()
if [[ -n "${NANOSIM_MODEL}" ]]; then
  nanosim_args=(--flags nanosim "--model ${NANOSIM_MODEL}")
fi

samovar prepare \
  --input_dir "${READS}" \
  --output_dir "${OUT}" \
  --cores "${CORES}" \
  --max-genomes "${MAX_GENOMES}" \
  --max-genome-mb "${MAX_GENOME_MB}" \
  "${skip_args[@]}" \
  --N 1 \
  --N_reads "${N_READS}" \
  --table_reads_generator concat \
  --reads_generator nanosim \
  "${nanosim_args[@]}" \
  --qc chopper \
  --reprofiler ensemble \
  --metax "metax ${METAX_DB}" \
  --kaiju "kaiju ${KAIJU_DB}" \
  --kraken2 "kraken2 ${KRAKEN2_DB}" \
  --length "length ." \
  --kmer "kmer_encoder . -k 4"

samovar exec --output_dir "${OUT}"
```

Длина — отдельная фича `feat_length_*_length`: сумма оснований R1 и R2. Это не поле length из отчёта Kraken.

## 3. Объединить и заново обучить модель

`samovar merge` требует одинаковый набор аннотаторов. Из Illumina уходят те, которых нет в ONT (Centrifuge, anvi'o, GC). Из ONT уходит длина, если её не было в Illumina. Остаются MetaX, Kaiju, Kraken 2 и k-mer k=4.

Скрипт не трогает исходные запуски: делает hardlink-копии, переписывает в копиях конфиг и таблицы, потом мержит `regenerated` и стартует с `reprofile`.

```bash
#!/usr/bin/env bash
set -euo pipefail

# --- пути этого кластера ---
SAMOVAR=/path/to/samovar
ILLUMINA_OUT=/path/to/work/samovar_illumina
ONT_OUT=/path/to/work/samovar_ont
WORK=/path/to/work/samovar_merge_stage
MERGED=/path/to/work/samovar_illumina_ont

export PATH="${SAMOVAR}/bin:${PATH}"
export PYTHONPATH="${SAMOVAR}/src${PYTHONPATH:+:${PYTHONPATH}}"

rm -rf "${WORK}"
mkdir -p "${WORK}"
cp -al "${ILLUMINA_OUT}" "${WORK}/illumina"
cp -al "${ONT_OUT}" "${WORK}/ont"

python3 - "${WORK}/illumina" "${WORK}/ont" <<'PY'
import sys
from pathlib import Path

import pandas as pd
import yaml
from omegaconf import OmegaConf

from samovar.add_annotator import existing_run_names


def rewrite_text(path: Path, text: str) -> None:
    path.unlink()
    path.write_text(text, encoding="utf-8")


def drop_run_config(path: Path, keep: set) -> None:
    if not path.is_file():
        return
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    rows = []
    for row in data.get("run_config") or []:
        name = str((row or {}).get("run_name") or "").strip()
        if name in keep:
            rows.append(row)
    data["run_config"] = rows
    rewrite_text(path, yaml.safe_dump(data, sort_keys=False))


def drop_hydra_cmds(root: Path, keep: set) -> None:
    path = root / ".hydra" / "config.yaml"
    if not path.is_file():
        return
    snap = OmegaConf.load(path)
    args = snap.get("prepare")
    args = None if args is None else args.get("args")
    if args is None:
        return
    for key in list(args.keys()):
        if str(key).startswith("cmd_") and str(key)[4:] not in keep:
            del args[key]
    path.unlink()
    OmegaConf.save(snap, path)


def column_kept(name: str, keep: set) -> bool:
    text = str(name)
    low = text.lower()
    for run in keep:
        if low.startswith(f"taxid_{run.lower()}_") or low.startswith(f"feat_{run.lower()}_"):
            return True
    if low.startswith("taxid_") or low.startswith("feat_"):
        return False
    return True


def filter_annotations(root: Path, keep: set) -> None:
    for folder in ("initial_annotations", "regenerated_annotations"):
        ann = root / folder
        if not ann.is_dir():
            continue
        for csv in ann.glob("*.csv"):
            frame = pd.read_csv(csv)
            cols = [c for c in frame.columns if column_kept(c, keep)]
            if cols == list(frame.columns):
                continue
            csv.unlink()
            frame.loc[:, cols].to_csv(csv, index=False)


illumina, ont = Path(sys.argv[1]), Path(sys.argv[2])
keep = set(existing_run_names(illumina)) & set(existing_run_names(ont))
if not keep:
    raise SystemExit("no shared annotators")
print("shared annotators:", " ".join(sorted(keep)))
for root in (illumina, ont):
    drop_run_config(root / ".log" / "configs" / "config_init.yaml", keep)
    drop_run_config(root / ".log" / "configs" / "config_reannotate.yaml", keep)
    drop_hydra_cmds(root, keep)
    filter_annotations(root, keep)
PY

samovar merge --mode regenerated --force \
  --output_dir "${MERGED}" \
  "${WORK}/illumina" "${WORK}/ont"

samovar exec --output_dir "${MERGED}" --startpoint reprofile
```

После `exec` модель лежит в `${MERGED}/reprofiled_annotations/`.
