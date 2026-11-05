#!/usr/bin/env python3
"""Read-to-initial-MAG assignment with auditable alignment candidates."""
from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import sqlite3
import subprocess
from contextlib import contextmanager
from pathlib import Path

from samovar.MAG.homology import iter_mag_fastas
from samovar.seqio import iter_fastq_records, iter_seqio_fasta, strip_mate_suffix


def _fastq_records(path):
    if not path:
        return
    if not Path(path).is_file():
        raise FileNotFoundError(path)
    for record in iter_fastq_records(path):
        yield record[0].lstrip('@').split()[0], record[1].strip().upper()


def _tax_by_mag(path):
    from samovar.assembly_profiling import load_mag_taxonomy
    if not path:
        return {}
    if not Path(path).is_file():
        raise FileNotFoundError(path)
    return load_mag_taxonomy(Path(path))


def _classified_mags(path):
    if not path:
        return set()
    result = set()
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            lineage = str(row.get('lineage') or '').strip()
            taxid = str(row.get('taxid') or row.get('taxID') or '0')
            if lineage.lower().startswith(('unclassified', 'n/a')):
                continue
            if taxid != '0' or lineage.startswith('d__'):
                result.add(row['mag_id'])
    return result


def contained_mag(read, mags):
    """Exact-match baseline: an ambiguous or absent match has no MAG."""
    if len(read) < 8:
        return ''
    reverse = read.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
    hits = {name for name, sequence in mags if read in sequence or reverse in sequence}
    return next(iter(hits)) if len(hits) == 1 else ''


def map_reads(r1, r2, mag_dir, *, fallback=False, dest_for_fallback=None):
    """Small exact-match baseline; never invent an assignment for a missed read."""
    if fallback:
        raise ValueError('Identity fallback is not supported for MAG mapping')
    mags = [(f.stem, str(rec.seq).upper()) for f in iter_mag_fastas(Path(mag_dir))
            for rec in iter_seqio_fasta(f)]
    return [(name, seq, mag) for path in (r1, r2) for name, seq in _fastq_records(path)
            if (mag := contained_mag(seq, mags))]


@contextmanager
def sam_lines(path):
    from samovar.assembly_profiling import _alignment_needs_samtools, which_tool

    if not _alignment_needs_samtools(Path(path)):
        with Path(path).open() as handle:
            yield handle
        return
    with subprocess.Popen([which_tool('samtools'), 'view', str(path)], stdout=subprocess.PIPE,
                          text=True) as proc:
        try:
            yield proc.stdout
        finally:
            if proc.stdout:
                proc.stdout.close()
        code = proc.wait()
        if code:
            raise subprocess.CalledProcessError(code, proc.args)


def marker_intervals(path):
    if not path:
        return {}
    result = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            result.setdefault((row['mag_id'], row['contig']), []).append(
                (int(row['start']), int(row['end'])))
    return result


def alignment_candidate(line, markers):
    p = line.rstrip('\n').split('\t')
    if len(p) < 11:
        raise ValueError('Malformed SAM alignment')
    flag = int(p[1])
    if flag & 4 or p[2] == '*':
        return None
    if '~' not in p[2]:
        raise ValueError(f'Alignment reference lacks MAG namespace: {p[2]}')
    mag, contig = p[2].split('~', 1)
    tags = {v.split(':', 2)[0]: v.split(':', 2)[2] for v in p[11:] if v.count(':') >= 2}
    ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', p[5])]
    if not ops or ''.join(f'{n}{op}' for n, op in ops) != p[5]:
        raise ValueError('Malformed mapped CIGAR')
    start = int(p[3]) - 1
    pos = start
    blocks = []
    for n, op in ops:
        if op in 'M=X':
            blocks.append((pos, pos + n))
        if op in 'MDN=X':
            pos += n
    covered = sum(n for n, op in ops if op in 'MI=X')
    if 'NM' not in tags and any(op == 'M' for _, op in ops):
        raise ValueError('MAG mapping requires NM tags (minimap2 SAM output)')
    matches = max(0, sum(n for n, op in ops if op in 'MID=X') - int(tags.get('NM', 0)))
    # Union marker intervals prevents counting overlapping marker annotations twice.
    overlap = []
    for left, right in blocks:
        for a, b in markers.get((mag, contig), []):
            if max(left, a) < min(right, b):
                overlap.append((max(left, a), min(right, b)))
    marker_bases = 0
    end = -1
    for a, b in sorted(overlap):
        marker_bases += max(0, b - max(a, end))
        end = max(end, b)
    if 'AS' not in tags:
        raise ValueError('MAG mapping requires AS tags')
    mate = 1 if flag & 64 else (2 if flag & 128 else (int(p[0][-1]) if p[0].endswith(('/1', '/2')) else 0))
    return (strip_mate_suffix(p[0]), mate, mag, contig, start, pos, flag,
            int(tags['AS']), matches, covered, marker_bases, int(p[4]))


def candidate_vector(ranks):
    """Share of summed alignment score for each classified MAG, sorted by id."""
    if not ranks:
        return ''
    total = sum(rank[0] for rank in ranks.values())
    names = sorted(ranks)
    if total <= 0:
        share = 1.0 / len(names)
        return ','.join(f'{mag}={share:.6f}' for mag in names)
    return ','.join(f'{mag}={ranks[mag][0] / total:.6f}' for mag in names)


def assign_alignments(r1, r2, bam, taxonomy, dest, *, markers='', sample='', settings=None, keep_vector=False):
    """Persist every reported candidate and emit one resolved row per read id.

    Within each MAG use the best alignment per mate: supplementary/secondary
    records cannot inflate the score by recounting the same query bases.
    """
    out = Path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)
    database = out.with_suffix(out.suffix + '.candidates.sqlite')
    pending = database.with_suffix('.sqlite.pending')
    pending.unlink(missing_ok=True)
    tax = _tax_by_mag(taxonomy)
    classified = _classified_mags(taxonomy)
    intervals = marker_intervals(markers)
    counts = dict(total=0, unmapped=0, multi_hit=0, unresolved=0, unclassified=0)
    try:
        with sqlite3.connect(pending) as db:
            db.executescript('''
                CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT);
                CREATE TABLE reads (seq TEXT PRIMARY KEY);
                CREATE TABLE candidates (
                  seq TEXT, mate INTEGER, mag_id TEXT, contig TEXT,
                  start INTEGER, end INTEGER, flag INTEGER, score INTEGER,
                  matches INTEGER, covered INTEGER, marker_bases INTEGER, mapq INTEGER,
                  selected INTEGER DEFAULT 0);
                CREATE INDEX candidate_reads ON candidates(seq,mag_id,mate);
                CREATE TABLE assignments (seq TEXT PRIMARY KEY,taxid TEXT,mag_id TEXT,status TEXT);
            ''')
            db.executemany('INSERT INTO metadata VALUES (?,?)', [('schema_version','1'), ('sample',sample),
                ('alignment_settings', json.dumps(settings or {}, sort_keys=True))])
            for path in (r1, r2):
                db.executemany('INSERT OR IGNORE INTO reads VALUES (?)',
                               ((strip_mate_suffix(name),) for name, _ in _fastq_records(path)))
            if bam:
                with sam_lines(bam) as lines:
                    for line in lines:
                        if not line.strip() or line.startswith('@'):
                            continue
                        hit = alignment_candidate(line, intervals)
                        if hit:
                            if hit[2] not in tax:
                                raise ValueError(f'Alignment MAG missing from taxonomy: {hit[2]}')
                            if not db.execute('SELECT 1 FROM reads WHERE seq=?', (hit[0],)).fetchone():
                                raise ValueError(f'Alignment query missing from FASTQ: {hit[0]}')
                            db.execute('INSERT INTO candidates VALUES (?,?,?,?,?,?,?,?,?,?,?,?,0)', hit)
            tmp_out = out.with_suffix(out.suffix + '.pending')
            with tmp_out.open('w') as handle:
                handle.write('seq\ttaxID\tMAG_ID' + ('\tMAG_vector\n' if keep_vector else '\n'))
                for (seq,) in db.execute('SELECT seq FROM reads ORDER BY seq'):
                    best = {}
                    reported_mags = set()
                    for row in db.execute('SELECT mate,mag_id,score,matches,covered,marker_bases,mapq '
                                          'FROM candidates WHERE seq=?', (seq,)):
                        reported_mags.add(row[1])
                        if row[1] not in classified:
                            continue
                        key = (row[1], row[0])
                        best[key] = max(best.get(key, row[2:]), row[2:])
                    ranks = {}
                    for (mag, _), rank in best.items():
                        ranks[mag] = tuple(a+b for a,b in zip(ranks.get(mag, (0,)*5), rank))
                    winners = sorted(mag for mag, rank in ranks.items() if rank == max(ranks.values()))
                    mag = winners[0] if len(winners) == 1 else ''
                    tids = {tax[m] for m in winners}
                    tid = next(iter(tids)) if len(tids) == 1 else '0'
                    status = 'mapped' if mag else ('unresolved' if winners else ('unclassified' if reported_mags else 'unmapped'))
                    counts['total'] += 1
                    counts['multi_hit'] += len(reported_mags) > 1
                    counts['unresolved'] += status == 'unresolved'
                    counts['unmapped'] += status == 'unmapped'
                    counts['unclassified'] += status == 'unclassified'
                    db.execute('INSERT INTO assignments VALUES (?,?,?,?)', (seq,tid,mag,status))
                    if mag:
                        db.execute('UPDATE candidates SET selected=1 WHERE seq=? AND mag_id=?', (seq,mag))
                    vector = '\t' + candidate_vector(ranks) if keep_vector else ''
                    handle.write(f'{seq}\t{tid}\t{mag}{vector}\n')
            db.commit()
        pending.replace(database)
        tmp_out.replace(out)
    except BaseException:
        pending.unlink(missing_ok=True)
        out.with_suffix(out.suffix + '.pending').unlink(missing_ok=True)
        raise
    report = {**counts, **{key + '_fraction': counts[key] / max(1, counts['total'])
                          for key in ('unmapped', 'multi_hit', 'unresolved', 'unclassified')}}
    out.with_suffix(out.suffix + '.mapping.json').write_text(json.dumps(report, indent=2) + '\n')
    return out


def write_mapping_table(r1, r2, mag_dir, taxonomy, dest, *, feature='MAG_ID', threads=1, markers='', aligner='minimap2', keep_vector=False):
    from samovar.assembly_profiling import run_minimap2, _fasta_has_records

    if feature != 'MAG_ID':
        raise ValueError('Legacy MAP_MAG_ID requires rerunning initial and downstream annotation with MAG_ID')
    bam = ''
    flags = []
    if aligner != 'exact':
        from samovar.main_config import imported_flags_for_names, iter_tools
        from samovar.paths import load_config
        flags = shlex.split(imported_flags_for_names(iter_tools(load_config()), 'minimap2'))
        flags.append('--secondary=yes')
    if any(_fasta_has_records(p) for p in iter_mag_fastas(Path(mag_dir))):
        bam = str(Path(dest).with_suffix('.mapping.bam'))
        if aligner == 'exact':
            # Contract baseline only: real exact matches, never a fabricated hit.
            bam = str(Path(dest).with_suffix('.mapping.sam'))
            references = [(f.stem, rec.id, str(rec.seq).upper())
                          for f in iter_mag_fastas(Path(mag_dir)) for rec in iter_seqio_fasta(f)]
            with Path(bam).open('w') as handle:
                for path in (r1, r2):
                    for name, seq in _fastq_records(path):
                        reverse = seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
                        for mag, contig, reference in references:
                            for flag, query in ((0, seq), (16, reverse)):
                                start = reference.find(query)
                                if start < 0:
                                    continue
                                handle.write(f'{name}\t{flag}\t{mag}~{contig}\t{start+1}\t0\t{len(seq)}M\t*\t0\t0\t{query}\t*\tAS:i:{len(seq)}\tNM:i:0\n')
        else:
            run_minimap2(r1, r2, mag_dir, bam, threads=threads, extra=flags)
    return assign_alignments(r1, r2, bam, taxonomy, dest, markers=markers, sample=Path(dest).stem,
                             keep_vector=keep_vector,
                             settings={'aligner':aligner, 'secondary':True, 'threads':threads,
                                       'presets_by_library':{'short':'sr','ont':'map-ont'}, 'flags':flags,
                                       'candidate_limit':'aligner defaults plus imported flags',
                                       'selection':'AS,matches,covered,marker_bases,MAPQ; best per mate/MAG',
                                       'candidate_vector': 'AS share across classified MAGs' if keep_vector else ''})


def validate_mapping_resume(output_dir):
    """Reject legacy ML feature schemas rather than silently renaming them."""
    root = Path(output_dir)
    for folder in ('initial_annotations', 'regenerated_annotations'):
        parent = root / folder
        if not parent.is_dir():
            continue
        for path in sorted(parent.glob('*.csv')):
            with path.open() as handle:
                header = next(csv.reader(handle), [])
            if any('MAP_MAG_ID' in name for name in header):
                raise ValueError('Legacy MAP_MAG_ID annotation: prepare a fresh run and rerun initial annotation and downstream stages; columns cannot be renamed safely')


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True)
    parser.add_argument('-I', dest='r2', default='')
    parser.add_argument('-r', required=True)
    parser.add_argument('-x', required=True)
    parser.add_argument('-o', required=True)
    parser.add_argument('-t', type=int, default=1)
    parser.add_argument('--markers', default='')
    parser.add_argument('--aligner', choices=('exact','minimap2'), default='exact',
                        help='Exact-match contract baseline; production uses minimap2')
    parser.add_argument('--keep-candidate-vector', action='store_true')
    args = parser.parse_args(argv)
    write_mapping_table(args.i,args.r2,args.r,args.x,args.o,threads=args.t,markers=args.markers,
                        aligner=args.aligner, keep_vector=args.keep_candidate_vector)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
