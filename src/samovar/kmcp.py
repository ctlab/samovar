"""Native KMCP search/profile-binning adapter for the read annotation contract."""
from __future__ import annotations

import argparse
import csv
import json
import shlex
import subprocess
import sys
import tempfile
from itertools import chain, zip_longest
from pathlib import Path

from samovar.paths import command_argv, resolve_executable
from samovar.seqio import fastq_record_id, iter_fastq_records


def database_files(directory):
    root = Path(directory)
    index_info = next((root / 'index').rglob('__db.yml'), root / 'index' / '__db.yml')
    paths = [index_info, root / 'taxid.map',
             root / 'taxdump' / 'nodes.dmp', root / 'taxdump' / 'names.dmp']
    for path in paths:
        if not path.is_file() or not path.stat().st_size:
            raise FileNotFoundError(f'KMCP database requires {path}')
    mapping = {}
    for line in paths[1].read_text().splitlines():
        fields = line.split('\t')
        if len(fields) != 2 or not fields[0] or not fields[1].isdigit() or int(fields[1]) < 1 or fields[0] in mapping:
            raise ValueError('KMCP taxid.map requires unique reference IDs and positive NCBI taxids')
        mapping[fields[0]] = fields[1]
    if not mapping:
        raise ValueError('Empty KMCP taxid.map')
    return root / 'index', paths[1], root / 'taxdump'


def input_ids(r1, r2=None):
    """Validate both files before native execution; pair IDs use SAMOVAR semantics."""
    mates = iter_fastq_records(r2) if r2 else iter(())
    first = next(mates, None)
    paired = first is not None
    records = zip_longest(iter_fastq_records(r1), chain([first], mates)) if paired else (
        (record, None) for record in iter_fastq_records(r1))
    ids, seen = [], set()
    for left, right in records:
        if left is None or (paired and right is None):
            raise ValueError('KMCP paired FASTQ record counts differ')
        for record in (left, right) if paired else (left,):
            seq, qual = record[1].strip(), record[3].rstrip('\r\n')
            if not seq or len(seq) != len(qual) or set(seq.upper()) - set('ACGTRYSWKMBDHVN'):
                raise ValueError('Malformed FASTQ sequence/quality for KMCP')
        name = fastq_record_id(left[0])
        if paired and fastq_record_id(right[0]) != name:
            raise ValueError('KMCP paired FASTQ identifiers differ')
        if name in seen:
            raise ValueError(f'Duplicate KMCP input read identifier: {name}')
        seen.add(name)
        ids.append(name)
    return ids, paired


def read_binning(path):
    """Consume CAMI assignments only; a search target is never treated as a taxid."""
    assignments = {}
    header = False
    for line in Path(path).read_text().splitlines():
        if line == '@@SEQUENCEID\tTAXID':
            header = True
            continue
        if not line or line.startswith(('#', '@')):
            continue
        fields = line.split('\t')
        if not header or len(fields) != 2 or not fields[1].isdigit():
            raise ValueError('Malformed KMCP CAMI binning result')
        name, taxid = fastq_record_id(fields[0]), str(int(fields[1]))
        if name in assignments and assignments[name] != taxid:
            raise ValueError(f'Conflicting KMCP assignments for {name}')
        assignments[name] = taxid
    if not header:
        raise ValueError('Missing KMCP CAMI binning header')
    return assignments


def native_flags(text, reserved):
    flags = shlex.split(text or '')
    for token in flags:
        if token.split('=')[0] in reserved or (token.startswith('-') and token[:2] in reserved):
            raise ValueError(f'KMCP adapter owns input/output/database/thread option: {token}')
    return flags


def annotate(r1, r2, database, output, threads=1, executable='kmcp', search_flags='', profile_flags=''):
    index, mapping, taxdump = database_files(database)
    cmd = command_argv(resolve_executable(executable, tool_key='kmcp'))
    if not cmd or not Path(cmd[0]).is_file():
        raise FileNotFoundError('KMCP is not installed; use ./install.sh KMCP')
    ids, paired = input_ids(r1, r2)
    common = {'-o', '--out-file', '-j', '--threads', '-i', '--infile-list'}
    search_extra = native_flags(search_flags, common | {'-d', '--db-dir', '-1', '-2', '--read1', '--read2',
                                                      '-N', '--name-map', '-D', '--default-name-map', '-g', '--query-whole-file'})
    profile_extra = native_flags(profile_flags, common | {'-T', '--taxid-map', '-X', '--taxdump', '-B', '--binning-result'})
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.parent / '.tmp' / 'kmcp'
    tmp.mkdir(parents=True, exist_ok=True)
    commands = []

    def run(argv):
        commands.append(argv)
        print(shlex.join(argv), file=sys.stderr)
        result = subprocess.run(argv, capture_output=True, text=True)
        print(result.stdout, end='', file=sys.stderr)
        print(result.stderr, end='', file=sys.stderr)
        result.check_returncode()
        return result

    assignments = {}
    with tempfile.TemporaryDirectory(dir=tmp, prefix='annotation-') as work:
        work = Path(work)
        version = run([*cmd, 'version']).stdout.strip()
        if ids:
            search, profile, binning = work / 'search.tsv', work / 'profile.tsv', work / 'reads.binning'
            reads = ['-1', str(r1), '-2', str(r2)] if paired else [str(r1)]
            run([*cmd, 'search', '-d', str(index), '-j', str(max(1, int(threads))), '-o', str(search), *search_extra, *reads])
            if not search.is_file():
                raise RuntimeError('KMCP search succeeded without its output')
            result = run([*cmd, 'profile', '-T', str(mapping), '-X', str(taxdump), '-B', str(binning),
                          '-j', str(max(1, int(threads))), '-o', str(profile), *profile_extra, str(search)])
            if binning.is_file():
                assignments = read_binning(binning)
            elif 'no target found, please adjust parameters' not in result.stderr:
                raise RuntimeError('KMCP profile succeeded without binning or an explicit no-target result')
        unknown = assignments.keys() - set(ids)
        if unknown:
            raise ValueError(f'KMCP returned unknown read identifiers: {sorted(unknown)[:5]}')
        pending = work / 'annotation.tsv'
        with pending.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['seq', 'taxID'])
            writer.writerows((name, assignments.get(name, '0')) for name in ids)
        pending.replace(out)
    Path(str(out) + '.kmcp.json').write_text(json.dumps(dict(
        version=version, commands=commands, database=str(Path(database).resolve()),
        input_records=len(ids), assigned_records=sum(t != '0' for t in assignments.values()),
        paired=paired, unclassified='0'), indent=2) + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', required=True)
    parser.add_argument('-I', default='')
    parser.add_argument('-d', required=True)
    parser.add_argument('-o', required=True)
    parser.add_argument('-t', type=int, default=1)
    parser.add_argument('--exec-path', default='kmcp')
    parser.add_argument('--search-flags', default='')
    parser.add_argument('--profile-flags', default='')
    args = parser.parse_args(argv)
    annotate(args.i, args.I, args.d, args.o, args.t, args.exec_path, args.search_flags, args.profile_flags)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
