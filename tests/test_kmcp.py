"""KMCP's native adapter boundaries; native ecosystem checks are separate."""
import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from samovar.kmcp import annotate, database_files, input_ids, read_binning


def database(tmp_path):
    for name in ('index/R001/__db.yml', 'taxid.map', 'taxdump/nodes.dmp', 'taxdump/names.dmp'):
        p = tmp_path / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('ref_562\t562\n' if name == 'taxid.map' else 'fixture')
    return tmp_path


def reads(path, names=('r/1', 'r/10')):
    path.write_text(''.join(f'@{name} description\nACGT\n+\nIIII\n' for name in names))
    return path


def test_kmcp_input_identity_and_pairs(tmp_path):
    r1 = reads(tmp_path / 'r1.fq', ('r/1', 'long/10'))
    r2 = reads(tmp_path / 'r2.fq', ('r/2', 'long/10'))
    assert input_ids(r1, r2) == (['r', 'long/10'], True)
    r2.write_text('')
    assert input_ids(r1, r2) == (['r', 'long/10'], False)
    reads(r2, ('other',))
    with pytest.raises(ValueError, match='identifiers differ'):
        input_ids(r1, r2)
    reads(r2, ('r/2',))
    with pytest.raises(ValueError, match='counts differ'):
        input_ids(r1, r2)
    r1.write_text('@r\nACGT\n+\nI\n')
    with pytest.raises(ValueError, match='sequence/quality'):
        input_ids(r1)
    r1.write_text('@r\nACGT\n+\n')
    with pytest.raises(ValueError, match='Truncated'):
        input_ids(r1)
    reads(r1, ('r/1', 'r/2'))
    with pytest.raises(ValueError, match='Duplicate'):
        input_ids(r1)


@pytest.mark.parametrize('missing', ['index/R001/__db.yml', 'taxid.map', 'taxdump/nodes.dmp', 'taxdump/names.dmp'])
def test_kmcp_missing_database_fails_even_with_empty_reads(tmp_path, missing):
    db = database(tmp_path / 'db')
    (db / missing).unlink()
    r1 = tmp_path / 'empty.fq'; r1.touch()
    with pytest.raises(FileNotFoundError, match=Path(missing).name):
        annotate(r1, None, db, tmp_path / 'out')


def test_kmcp_binning_is_required_and_strict(tmp_path):
    path = tmp_path / 'binning.tsv'
    path.write_text('@Version:0.10.0\n@@SEQUENCEID\tTAXID\nr/1\t562\nlong/10\t4932\n')
    assert read_binning(path) == {'r': '562', 'long/10': '4932'}
    path.write_text('@@SEQUENCEID\tTAXID\nr/1\t562\nr/2\t4932\n')
    with pytest.raises(ValueError, match='Conflicting'):
        read_binning(path)
    path.write_text('r\t126\t10\tref_562\n')
    with pytest.raises(ValueError, match='CAMI'):
        read_binning(path)
    db = database(tmp_path / 'db'); (db / 'taxid.map').write_text('ref_562\tGTDB_SGB1\n')
    with pytest.raises(ValueError, match='positive NCBI'):
        database_files(db)


def fake_native(monkeypatch, tmp_path, outcome='assigned'):
    binary = tmp_path / 'program with spaces' / 'kmcp'
    binary.parent.mkdir(); binary.touch()
    monkeypatch.setattr('samovar.kmcp.resolve_executable', lambda *a, **k: str(binary))
    calls = []
    def run(argv, **kwargs):
        calls.append(argv)
        stdout, stderr, rc = '', '', 0
        if 'version' in argv:
            stdout = 'kmcp v0.9.5\n'
        elif 'search' in argv:
            Path(argv[argv.index('-o')+1]).write_text('native search fixture')
        elif 'profile' in argv:
            if outcome == 'assigned':
                Path(argv[argv.index('-B')+1]).write_text('@@SEQUENCEID\tTAXID\nr/1\t562\n')
            elif outcome == 'nohit':
                stderr = '[WARN] no target found, please adjust parameters\n'
            elif outcome == 'failed':
                rc, stderr = 2, 'technical failure'
        return subprocess.CompletedProcess(argv, rc, stdout, stderr)
    monkeypatch.setattr('samovar.kmcp.subprocess.run', run)
    return calls


def test_kmcp_annotation_uses_binning_and_preserves_unclassified(tmp_path, monkeypatch):
    calls = fake_native(monkeypatch, tmp_path)
    db = database(tmp_path / 'db with spaces')
    r1, r2 = reads(tmp_path / 'reads 1.fq'), reads(tmp_path / 'reads 2.fq', ('r/2', 'r/10'))
    output = tmp_path / 'output dir' / 'annotation.tsv'
    annotate(r1, r2, db, output, 3, profile_flags='-m 1', search_flags='-w')
    frame = pd.read_csv(output, sep='\t')
    assert frame.to_dict('list') == {'seq': ['r', 'r/10'], 'taxID': [562, 0]}
    assert calls[1][calls[1].index('-j')+1] == '3'
    assert calls[1][calls[1].index('-1')+1] == str(r1) and '-2' in calls[1]
    assert '-T' in calls[2] and '-X' in calls[2] and '-B' in calls[2]
    assert not list((output.parent / '.tmp/kmcp').iterdir())
    assert json.loads(Path(str(output)+'.kmcp.json').read_text())['assigned_records'] == 1


@pytest.mark.parametrize('outcome', ['nohit', 'missing', 'failed'])
def test_kmcp_nohit_and_technical_errors_are_distinct(tmp_path, monkeypatch, outcome):
    fake_native(monkeypatch, tmp_path, outcome)
    db = database(tmp_path / 'db'); r1 = reads(tmp_path / 'reads.fq'); output = tmp_path / 'out'
    output.write_text('previous result')
    if outcome == 'nohit':
        annotate(r1, None, db, output)
        assert pd.read_csv(output, sep='\t').taxID.eq(0).all()
    else:
        error = subprocess.CalledProcessError if outcome == 'failed' else RuntimeError
        with pytest.raises(error):
            annotate(r1, None, db, output)
        assert output.read_text() == 'previous result'


def test_kmcp_empty_input_and_reserved_flags(tmp_path, monkeypatch):
    calls = fake_native(monkeypatch, tmp_path)
    db = database(tmp_path / 'db'); r1 = tmp_path / 'reads.fq'; r1.touch(); output = tmp_path / 'out'
    annotate(r1, None, db, output)
    assert pd.read_csv(output, sep='\t').empty and len(calls) == 1
    with pytest.raises(ValueError, match='owns'):
        annotate(r1, None, db, output, profile_flags='-Tother.map')


@pytest.mark.parametrize("genome_count", [1, 2])
def test_kmcp_builder_uses_shared_references_and_native_resources(tmp_path, monkeypatch, genome_count):
    from samovar.build_database import build_database_kmcp
    genome = tmp_path / '562.fna'; genome.write_text('>chromosome\nACGTACGT\n')
    calls = []
    def run(argv):
        calls.append(argv)
        if 'compute' in argv:
            source = Path(argv[argv.index('-I')+1]); assert (source/'ref_562_0.fna').read_text() == genome.read_text()
            assert len(list(source.glob('*.fna'))) == genome_count
        if 'index' in argv:
            index = Path(argv[argv.index('-O')+1])/'R001';index.mkdir(parents=True)
            (index/'__db.yml').write_text('native');(index/'block.uniki').write_bytes(b'index')
    def taxdump(db_path, taxonomy_path=None):
        p = Path(db_path);p.mkdir()
        for name in ('nodes.dmp','names.dmp'):(p/name).write_text('fixture')
    monkeypatch.setattr('samovar.build_database.run_command', run)
    monkeypatch.setattr('samovar.build_database.get_taxonomy_db', taxdump)
    dest = tmp_path/'db with spaces'
    genomes = {str(genome): '562'}
    if genome_count == 2:
        second = tmp_path / '562_other.fna'; second.write_text(genome.read_text())
        genomes[str(second)] = '562'
    build_database_kmcp(genomes, dest, threads=3)
    assert (dest/'taxid.map').read_text() == ''.join(f'ref_562_{i}\t562\n' for i in range(genome_count))
    assert all(cmd[cmd.index('-j')+1]=='3' for cmd in calls)
    manifest=json.loads((dest/'kmcp.manifest.json').read_text())
    assert manifest['parameters']['chunks']==10 and len(manifest['taxdump_sha256']['nodes.dmp'])==64
    assert not list((dest/'.tmp').iterdir())
    with pytest.raises(FileExistsError):
        build_database_kmcp({str(genome):'562'}, dest)

def test_kmcp_workflow_preserves_sample_and_run_names_with_underscores(tmp_path):
    import yaml
    from samovar.paths import command_argv, resolve_executable, workflow_dir
    inputs = tmp_path / 'reads';inputs.mkdir()
    reads(inputs / 'sample_one.fastq', ('r1',))
    cfg = dict(r1_dir=str(inputs), r2_dir=str(inputs), output_dir=str(tmp_path/'reports'),
               annotation_dir=str(tmp_path/'annotations'), cores=1,
               run_config=[dict(run_name='kmcp_named', type='kmcp', cmd='kmcp', db_path=str(tmp_path), threads=1)])
    config = tmp_path/'config.yaml';config.write_text(yaml.safe_dump(cfg))
    result = subprocess.run([*command_argv(resolve_executable('snakemake',tool_key='snakemake')),
                             '-s',str(Path(workflow_dir())/'annotators/Snakefile'),
                             '--configfile',str(config),'--cores','1','--dry-run'],capture_output=True,text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'sample=sample_one, run_name=kmcp_named' in result.stdout


def test_kmcp_named_run_preserves_sample_through_ingestion(tmp_path):
    from samovar.combine_tables import ingest_annotator, finalize_annotations
    out = tmp_path / 'sample_one_kmcp_native.custom_kmcp.out'
    out.write_text('seq\ttaxID\nread/10\t562\n')
    ann = tmp_path / 'annotations'
    path = ingest_annotator(out, ann, tmp_path / 'ingested', sample='sample_one')
    assert path.name == 'sample_one.annotation.csv'
    finalize_annotations(ann)
    frame = pd.read_csv(ann / 'combined_annotation_table.csv')
    assert frame['sample'].tolist() == ['sample_one']
    assert frame['seq'].tolist() == ['read/10']
    assert frame['taxID_kmcp_0'].tolist() == [562]
