"""Phanta's native workflow adapter and typed species profile boundary."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from samovar.phanta import database_files, run
from samovar.profiles import load_profile


HEADER = 'name\ttaxonomy_id\tnew_est_reads\trel_taxon_abundance\n'
PROFILE = HEADER + 'Escherichia coli\t562\t80\t0.8\nnovel native species\t9000001\t20\t0.2\n'


def test_phanta_native_species_and_unresolved_ids(tmp_path):
    source = tmp_path / 'profile.phanta.tsv'
    source.write_text(PROFILE)
    table, info = load_profile(source)
    assert table.taxid.tolist() == ['562', 'unmapped:phanta:9000001']
    assert table.N_sample.tolist() == [0.8, 0.2]
    assert (info['unit'], info['basis'], info['rank'], info['normalization']) == ('fraction', 'organisms', 'species', 'classified')
    assert info['native_taxa']['unmapped:phanta:9000001']['name'] == 'novel native species'
    # A native ID colliding with an unrelated NCBI species is not silently trusted.
    source.write_text(PROFILE.replace('9000001', '4932'))
    assert load_profile(source)[0].taxid.tolist()[-1] == 'unmapped:phanta:4932'
    source.write_text(HEADER)
    assert load_profile(source)[0].empty


@pytest.mark.parametrize('text', [
    'not a profile\n', HEADER + 'x\tnative\t1\t1\n',
    HEADER + 'x\t562\t-1\t1\n', HEADER + 'x\t562\t1\tNaN\n',
    HEADER + 'x\t562\t1\t2\n', PROFILE + 'duplicate\t562\t1\t0\n',
])
def test_phanta_invalid_native_profiles_fail(tmp_path, text):
    path = tmp_path / 'native.tsv'
    path.write_text(text)
    with pytest.raises((ValueError, KeyError)):
        load_profile(path, {'format': 'phanta'})


@pytest.fixture
def native(tmp_path, monkeypatch):
    root = tmp_path / "native space's Ω"
    root.mkdir()
    source = root / 'upstream'
    source.mkdir()
    executable = source / 'run_phanta.py'
    executable.write_text('')
    (source / 'Snakefile').write_text('# native workflow')
    (source / 'config.yaml').write_text('read_length: 150\nclass_threads: 16\n')
    (source / 'pipeline_scripts').mkdir()
    (source / 'pipeline_scripts/scale_bracken.py').write_text(
        'import pandas as pd\ngenome_length_path = "database/library/species_genome_size.txt"\n'
        'table = pd.read_csv(genome_length_path, sep="\\t")\n'
        'assert table.species_level_taxa.tolist() == ["562"]\n')
    db = root / 'database'
    required = ['hash.k2d', 'opts.k2d', 'taxo.k2d', 'inspect.out', 'taxonomy/nodes.dmp',
                'taxonomy/names.dmp', 'library/species_genome_size.txt', 'database126mers.kmer_distrib']
    for name in required:
        path = db / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('native database file\n')
    (db / 'library/species_genome_size.txt').write_text('species_level_taxa\tmedian\n562\t4000000\n')
    fake = root / 'snakemake'
    fake.write_text(f'#!{sys.executable}\n' + f'''import os, sys, subprocess
from pathlib import Path
import yaml
cfg = yaml.safe_load(Path('config.yaml').read_text())
assert cfg['database'] == 'database' and cfg['read_length'] == 126
assert cfg['class_threads'] == int(sys.argv[sys.argv.index('--cores') + 1])
assert Path(cfg['pipeline_directory'], 'Snakefile').is_file()
inputs = Path('samples.tsv').read_text().strip().split('\\t')[1:]
assert all(Path(p).is_file() and ' ' not in p for p in inputs)
subprocess.run(['python', '-c', 'import sys; assert sys.executable == sys.argv[1]', {sys.executable!r}], check=True)
db = Path('database')
if (db/'fail').exists(): raise SystemExit(7)
target = Path(sys.argv[-1]); target.parent.mkdir(parents=True, exist_ok=True)
if str(target).endswith('.krak.report'):
    target.write_text('100\\t100\\t100\\t10\\t10\\tS\\t562\\tEscherichia coli\\n' if not (db/'empty').exists()
                      else '100\\t100\\t100\\t0\\t0\\tU\\t0\\tunclassified\\n')
else:
    subprocess.run(['python', 'phanta/pipeline_scripts/scale_bracken.py'], check=True)
    target.write_text({PROFILE!r})
''')
    fake.chmod(0o755)
    cfg_path = Path(os.environ['SAMOVAR_CONFIG'])
    cfg = json.loads(cfg_path.read_text())
    replaced = {'phanta', 'kraken2', 'bracken', 'snakemake'}
    cfg['tools'] = {k: v for k, v in cfg['tools'].items() if k.split(':')[0] not in replaced}
    cfg['tool_envs'] = {k: v for k, v in cfg.get('tool_envs', {}).items() if k not in replaced}
    for name in ('phanta', 'kraken2', 'bracken', 'snakemake'):
        cfg['tools'][name] = {'type': 'profiler' if name == 'phanta' else 'workflow',
                              'exec': {'path': str(executable if name == 'phanta' else fake)}}
    cfg_path.write_text(json.dumps(cfg))
    foreign = root / 'foreign'
    foreign.mkdir()
    for name in ('python', 'python3', 'kraken2', 'bracken', 'snakemake'):
        path = foreign / name
        path.write_text('#!/bin/sh\nexit 99\n')
        path.chmod(0o755)
    monkeypatch.setenv('PATH', str(foreign) + os.pathsep + os.environ['PATH'])
    reads = root / 'R1.fastq'
    reads.write_text('@read\nACGT\n+\nIIII\n')
    return root, executable, db, reads


@pytest.mark.parametrize('paired', [False, True])
def test_phanta_production_profile_adapter(native, paired):
    from samovar.annotators_wrapper import ProfileAnnotator
    root, executable, db, reads = native
    tool = ProfileAnnotator(dict(type='phanta', cmd='phanta', run_name='native_profile', db_path=str(db), threads=2), {})
    output = tool.get_expected_outputs('patient_A', str(root))
    command = tool.get_snakemake_shell_cmd(str(reads), str(reads) if paired else '', output)
    subprocess.run(['bash', '-e', '-c', command], check=True)
    table, info = load_profile(output[0])
    assert table.N_patient_A.tolist() == [0.8, 0.2]
    assert info['identity'] == 'native_profile' and info['database'] == str(db)
    assert info['native_settings']['read_length'] == 126
    assert not list((root / '.tmp/phanta').iterdir())


@pytest.mark.parametrize('mode', ['empty_input', 'empty_result', 'failure', 'missing_distribution'])
def test_phanta_empty_results_and_technical_errors(native, mode):
    root, executable, db, reads = native
    out = root / 'native.tsv'
    if mode == 'empty_input':
        reads.write_text('')
    elif mode == 'empty_result':
        (db / 'empty').touch()
    elif mode == 'failure':
        (db / 'fail').touch()
    else:
        (db / 'database126mers.kmer_distrib').unlink()
    if mode in {'failure', 'missing_distribution'}:
        out.write_text('previous result')
        with pytest.raises((subprocess.CalledProcessError, FileNotFoundError)):
            run(str(reads), '', db, out)
        assert out.read_text() == 'previous result'
    else:
        run(str(reads), '', db, out)
        assert load_profile(out, {'format': 'phanta'})[0].empty
    scratch = root / '.tmp/phanta'
    assert not scratch.exists() or not list(scratch.iterdir())
