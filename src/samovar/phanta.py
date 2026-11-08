"""Phanta's native per-sample workflow behind the sample-profile contract."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

import yaml

from samovar.paths import command_argv, python_path, resolve_executable, _existing_tool
from samovar.profiles import load_profile, metadata
from samovar.seqio import fastq_has_reads, is_gzip_file


def database_files(database, read_length=None):
    db = Path(database).resolve()
    lengths = [int(p.name[len('database'):-len('mers.kmer_distrib')])
               for p in db.glob('database*mers.kmer_distrib')
               if p.name[len('database'):-len('mers.kmer_distrib')].isdigit()]
    length = read_length or (lengths[0] if len(lengths) == 1 else 150)
    required = ['hash.k2d', 'opts.k2d', 'taxo.k2d', 'inspect.out',
                'taxonomy/nodes.dmp', 'taxonomy/names.dmp',
                'library/species_genome_size.txt', f'database{length}mers.kmer_distrib']
    for filename in required:
        path = db / filename
        if not path.is_file() or not path.stat().st_size:
            raise FileNotFoundError(f'Phanta database requires {path}; select --read-length for a matching Bracken distribution')
    return db, length


def run(r1, r2, database, output, sample='sample', threads=1, executable='phanta', flags=''):
    options = argparse.ArgumentParser(add_help=False)
    options.add_argument('--read-length', type=int)
    options.add_argument('--confidence-threshold', type=float)
    options.add_argument('--filter-thresh', type=int)
    for domain in ('bacterial', 'viral', 'arc', 'euk'):
        options.add_argument('--cov-thresh-' + domain, type=float)
        options.add_argument('--minimizer-thresh-' + domain, type=int)
    selected = vars(options.parse_args(shlex.split(flags)))
    if threads < 1 or (selected['read_length'] is not None and selected['read_length'] < 1):
        raise ValueError('Phanta requires positive threads and read length')
    for key, value in selected.items():
        if value is not None and (value < 0 or (key.startswith(('cov_', 'confidence_')) and not 0 <= value <= 1)):
            raise ValueError(f'Invalid Phanta {key}: {value}')
    db, length = database_files(database, selected['read_length'])
    native = command_argv(resolve_executable(executable, tool_key='phanta', allow_python_script=True))
    source = Path(native[-1]).resolve().parent
    if not (source / 'Snakefile').is_file() or not (source / 'config.yaml').is_file():
        raise FileNotFoundError('Phanta source missing; use ./install.sh Phanta or import its run_phanta.py')
    commands = {name: command_argv(resolve_executable(name, tool_key=name)) for name in ('kraken2', 'bracken', 'snakemake')}
    if _existing_tool('phanta-kraken2'):
        commands['kraken2'] = command_argv(resolve_executable('phanta-kraken2', tool_key='phanta-kraken2'))
    commands['python'] = [python_path()]
    for name, command in commands.items():
        if not command or not Path(command[0]).is_file():
            raise FileNotFoundError(f'Phanta requires configured {name}; use install.sh')
    out = Path(output).resolve()
    scratch = out.parent / '.tmp' / 'phanta'
    scratch.mkdir(parents=True, exist_ok=True)
    info = metadata(dict(unit='fraction', basis='organisms', rank='species', normalization='classified',
                         native_taxonomy_path=str(db / 'taxonomy')))
    with tempfile.TemporaryDirectory(prefix='profile_', dir=scratch) as tmp:
        work = Path(tmp)
        # Upstream shell templates do not quote paths: all native paths are relative.
        native_view = work / 'phanta'
        native_view.mkdir()
        (native_view / 'Snakefile').symlink_to(source / 'Snakefile')
        shutil.copytree(source / 'pipeline_scripts', native_view / 'pipeline_scripts')
        scale = native_view / 'pipeline_scripts/scale_bracken.py'
        # Native code converts Bracken IDs to strings but lets pandas infer the
        # genome-size IDs as integers in an entirely numeric taxonomy.
        scale_code = scale.read_text()
        typed_ids = 'pd.read_csv(genome_length_path, sep="\\t", dtype={"species_level_taxa": str})'
        original = 'pd.read_csv(genome_length_path, sep="\\t")'
        if original in scale_code:
            scale.write_text(scale_code.replace(original, typed_ids))
            info['native_compatibility'] = ['read genome-size species IDs as strings']
        (work / 'database').symlink_to(db, target_is_directory=True)
        inputs = []
        for number, filename in enumerate([r1, r2] if r2 and fastq_has_reads(r2) else [r1]):
            name = f'reads{number}.fastq' + ('.gz' if is_gzip_file(filename) else '')
            (work / name).symlink_to(Path(filename).resolve())
            inputs.append(name)
        if len({is_gzip_file(p) for p in (r1, r2) if p}) > 1:
            raise ValueError('Phanta requires both mates to use the same compression')
        config = yaml.safe_load((source / 'config.yaml').read_text())
        config.update({k: v for k, v in selected.items() if v is not None})
        config.update(pipeline_directory='phanta', sample_file='samples.tsv', outdir='result',
                      database='database', class_threads=threads, read_length=length,
                      gzipped=is_gzip_file(r1), single_end_krak=False, delete_intermediate=False)
        (work / 'config.yaml').write_text(yaml.safe_dump(config))
        (work / 'samples.tsv').write_text('\t'.join(['sample', *inputs]) + '\n')
        bin_dir = work / 'bin'
        bin_dir.mkdir()
        for name in ('python', 'kraken2', 'bracken'):
            launcher = bin_dir / name
            launcher.write_text('#!/bin/sh\nexec ' + shlex.join(commands[name]) + ' "$@"\n')
            launcher.chmod(0o755)
        env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ.get('PATH', ''))
        target = Path('result/classification/sample.krak.report.filtered.bracken.scaled')
        if fastq_has_reads(r1):
            argv = [*commands['snakemake'], '-s', 'phanta/Snakefile', '--configfile', 'config.yaml',
                    '--cores', str(threads), '--jobs', '1']
            report = Path('result/classification/sample.krak.report')
            subprocess.run([*argv, str(report)], cwd=work, env=env, check=True)
            rows = [line.split('\t') for line in (work / report).read_text().splitlines()]
            if not rows or any(len(row) != 8 or not row[1].isdigit() or not row[6].isdigit() for row in rows):
                raise ValueError('Malformed native Phanta Kraken2 minimizer report')
            # Upstream's empty-species pandas.apply fails before its no-hit handler.
            # No species assignments cannot produce a species abundance estimate.
            if any(row[5].startswith('S') and int(row[1]) > 0 for row in rows):
                subprocess.run([*argv, str(target)], cwd=work, env=env, check=True)
                produced = work / target
            else:
                produced = work / 'empty.tsv'
                produced.write_text('name\ttaxonomy_id\tnew_est_reads\trel_taxon_abundance\n')
                info['native_empty_reason'] = 'no species assignments in native Kraken2 report'
            load_profile(produced, {**info, 'format': 'phanta', 'sample': sample})
        else:
            produced = work / 'empty.tsv'
            produced.write_text('name\ttaxonomy_id\tnew_est_reads\trel_taxon_abundance\n')
        info.update(native_settings=config, native_commands=commands)
        produced.replace(out)
        Path(str(out) + '.profile.json').write_text(json.dumps(info, indent=2) + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', required=True)
    parser.add_argument('-I', default='')
    parser.add_argument('-d', required=True)
    parser.add_argument('-o', required=True)
    parser.add_argument('-t', type=int, default=1)
    parser.add_argument('--sample', default='sample')
    parser.add_argument('--exec-path', default='phanta')
    parser.add_argument('--native-flags', default='')
    args = parser.parse_args(argv)
    run(args.i, args.I, args.d, args.o, args.sample, args.t, args.exec_path, args.native_flags)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
