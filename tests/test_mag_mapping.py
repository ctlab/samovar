"""Stable initial MAG mapping, actual evidence, and cohort catalog contracts."""
import json
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from samovar.MAG.mapping import assign_alignments, map_reads
from samovar.MAG.homology import (
    extract_gtdb_markers, finalize_mapping_catalog, validate_mapping_catalog,
)


def fastq(path, names):
    path.write_text(''.join(f'@{name}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n' for name in names))
    return str(path)


def alignment(name, mag, *, score=12, matches=12, mapq=0, flag=0, start=1):
    return f'{name}\t{flag}\t{mag}~c\t{start}\t{mapq}\t12M\t*\t0\t0\tACGTACGTACGT\tIIIIIIIIIIII\tAS:i:{score}\tNM:i:{12-matches}\n'


def resolve(tmp_path, sam, taxa=('562','9606'), markers=''):
    reads = fastq(tmp_path / 'reads.fastq', ['unique', 'tie', 'miss'])
    tax = tmp_path / 'tax.tsv'
    tax.write_text('mag_id\ttaxid\tlineage\nM1\t'+taxa[0]+'\t\nM2\t'+taxa[1]+'\t\n')
    bam = tmp_path / 'reads.sam'; bam.write_text(sam)
    out = tmp_path / 'mapped.tsv'
    assign_alignments(reads, '', str(bam), str(tax), str(out), markers=markers, sample='sample')
    return pd.read_table(out, keep_default_na=False).set_index('seq'), out


def test_candidate_vector_is_optional_and_sums_to_one(tmp_path):
    from samovar.MAG.mapping import assign_alignments
    from samovar.assembly_profiling import parse_slot_extra

    assert parse_slot_extra('--mag-mode map-initial --keep-candidate-vector')['keep_candidate_vector'] is True
    reads = fastq(tmp_path / 'reads.fastq', ['unique', 'tie', 'miss'])
    tax = tmp_path / 'tax.tsv'
    tax.write_text('mag_id\ttaxid\tlineage\nM1\t562\t\nM2\t9606\t\n')
    sam = tmp_path / 'reads.sam'
    sam.write_text(alignment('unique', 'M1', score=20) + alignment('unique', 'M2', score=10) +
                   alignment('tie', 'M1', score=8) + alignment('tie', 'M2', score=8))
    out = tmp_path / 'mapped.tsv'
    assign_alignments(reads, '', str(sam), str(tax), str(out), keep_vector=True)
    table = pd.read_table(out, keep_default_na=False).set_index('seq')
    assert table.loc['unique', 'MAG_ID'] == 'M1'
    assert table.loc['unique', 'MAG_vector'] == 'M1=0.666667,M2=0.333333'
    assert table.loc['tie', 'MAG_ID'] == ''
    assert table.loc['tie', 'MAG_vector'] == 'M1=0.500000,M2=0.500000'
    assert table.loc['miss', 'MAG_vector'] == ''
    plain = tmp_path / 'plain.tsv'
    assign_alignments(reads, '', str(sam), str(tax), str(plain))
    assert plain.read_text().splitlines()[0] == 'seq\ttaxID\tMAG_ID'


def test_covered_phage_maps_and_sparse_genome_does_not(tmp_path):
    """Regenerated reads map only where the initial MAG actually contains sequence.

    The toy identity MAG is the concatenated initial reads, not the source
    genomes. A short phage that is fully stored maps; new reads from the
    uncovered part of a long genome stay unmapped. That is why a toy
    map-initial run can leave most regenerated reads unmapped while mapping
    every initial read.
    """
    from samovar.MAG.mapping import write_mapping_table
    from samovar.MAG.tetra import read_fasta_sequence

    root = Path(__file__).resolve().parents[1]
    phage = read_fasta_sequence(root / 'data/test_genomes/meta/2886930.fna')
    coli = read_fasta_sequence(root / 'data/test_genomes/meta/562.fna')
    mags = tmp_path / 'mags'
    mags.mkdir()
    (mags / 'phage.fa').write_text(f'>phi\n{phage}\n')
    (mags / 'coli.fa').write_text(f'>slice\n{coli[:120]}\n')
    tax = tmp_path / 'tax.tsv'
    tax.write_text('mag_id\ttaxid\tlineage\nphage\t2886930\td__Viruses\ncoli\t562\td__Bacteria\n')
    records = [(f'phage{i}', phage[start:start + 80]) for i, start in enumerate(range(0, len(phage) - 80, 400))]
    records += [(f'coli{i}', coli[start:start + 80]) for i, start in enumerate(range(200, 4000, 400))]
    records.append(('coli_seen', coli[:80]))
    lines = []
    for name, sequence in records:
        lines.extend([f'@{name}', sequence, '+', 'I' * len(sequence)])
    (tmp_path / 'regen.fastq').write_text('\n'.join(lines) + '\n')
    out = tmp_path / 'mapped.tsv'
    write_mapping_table(str(tmp_path / 'regen.fastq'), '', str(mags), str(tax), str(out), aligner='exact')
    table = pd.read_table(out, keep_default_na=False).set_index('seq')
    phage_ids = [name for name, _ in records if name.startswith('phage')]
    missed = [name for name, _ in records if name.startswith('coli') and name != 'coli_seen']
    assert (table.loc[phage_ids, 'MAG_ID'] == 'phage').all()
    assert (table.loc[missed, 'MAG_ID'] == '').all()
    assert table.loc['coli_seen', 'MAG_ID'] == 'coli'
    report = json.loads(Path(str(out) + '.mapping.json').read_text())
    assert report['unmapped'] == len(missed)
    assert report['unmapped'] > 0


def test_conflicts_candidates_and_zero_rows(tmp_path):
    table, out = resolve(tmp_path, alignment('unique','M1',score=20) + alignment('unique','M2') +
                         alignment('tie','M2',flag=256) + alignment('tie','M1'))
    assert table.loc['unique','MAG_ID'] == 'M1'
    assert str(table.loc['unique','taxID']) == '562'
    assert table.loc['tie','MAG_ID'] == table.loc['miss','MAG_ID'] == ''
    assert str(table.loc['tie','taxID']) == str(table.loc['miss','taxID']) == '0'
    with sqlite3.connect(str(out)+'.candidates.sqlite') as db:
        assert db.execute('select count(*) from candidates').fetchone()[0] == 4
        assert db.execute('select value from metadata where key="sample"').fetchone()[0] == 'sample'
    report = json.loads(Path(str(out)+'.mapping.json').read_text())
    assert report['multi_hit'] == 2 and report['unresolved'] == 1 and report['unmapped'] == 1


def test_tied_mags_share_taxon(tmp_path):
    table, _ = resolve(tmp_path, alignment('tie','M2')+alignment('tie','M1'), taxa=('562','562'))
    assert str(table.loc['tie','taxID']) == '562' and table.loc['tie','MAG_ID'] == ''


def test_markers_break_tie_not_secondary_duplication(tmp_path):
    markers = tmp_path / 'markers.tsv'
    markers.write_text('mag_id\tcontig\tstart\tend\nM2\tc\t0\t12\nM2\tc\t0\t12\n')
    table, _ = resolve(tmp_path, alignment('tie','M1')*3 + alignment('tie','M2',flag=16), markers=str(markers))
    assert table.loc['tie','MAG_ID'] == 'M2'


def test_mate_secondary_and_supplementary_do_not_inflate_scores(tmp_path):
    table, _ = resolve(tmp_path, alignment('unique/1','M1',flag=64,score=20) +
                       alignment('unique/2','M1',flag=128,score=20) +
                       alignment('unique/1','M2',score=30)*3)
    assert table.loc['unique','MAG_ID'] == 'M1'


def test_no_fictional_identity_fallback(tmp_path):
    mags = tmp_path / 'mags';mags.mkdir()
    (mags/'M1.fa').write_text('>c\nAAAAAAAAAAAA\n')
    reads = fastq(tmp_path/'reads.fastq', ['miss'])
    assert map_reads(reads,'',mags) == []
    with pytest.raises(ValueError,match='fallback'):
        map_reads(reads,'',mags,fallback=True)


def test_empty_catalog_manifest_and_corruption(tmp_path):
    catalog = tmp_path / 'catalog.tsv'
    with pytest.raises(ValueError,match='incomplete'):
        validate_mapping_catalog(catalog)
    finalize_mapping_catalog(catalog, samples=['s2','s1'])
    assert validate_mapping_catalog(catalog)['samples'] == ['s1','s2']
    (tmp_path/'taxonomy.tsv').write_text('corrupt')
    with pytest.raises(ValueError,match='corrupt'):
        validate_mapping_catalog(catalog)


def test_gtdb_marker_join_uses_own_gene_coordinates(tmp_path):
    (tmp_path/'gtdbtk.bac120.markers_summary.tsv').write_text(
        'name\tlist_unique_genes\tlist_multiple_genes\tlist_multiple_unique_genes\nM1\tPF1\t\t\n')
    genes = tmp_path/'marker_genes/M1';genes.mkdir(parents=True)
    (genes/'M1_protein.faa').write_text('>c_1 # 2 # 12 # 1 # ID=1_1\nAAAA\n>c_2 # 20 # 30 # 1 # ID=1_2\nAAAA\n')
    (genes/'M1_pfam_tophit.tsv').write_text('Gene Id\tTop hits\nc_1\tPF1,1e-5,50\nc_2\tUNRELATED,1e-5,100\n')
    out = extract_gtdb_markers(tmp_path,tmp_path/'markers.tsv')
    lines = out.read_text().splitlines()
    assert len(lines) == 2 and lines[1] == 'M1\tc\t1\t12\tc_1\tPF1'


def test_mapping_error_does_not_publish_output(tmp_path):
    with pytest.raises(ValueError,match='MAG namespace'):
        resolve(tmp_path,alignment('unique','M1').replace('M1~c','c'))
    assert not (tmp_path/'mapped.tsv').exists()


def test_legacy_mapping_annotation_requires_fresh_initial(tmp_path):
    from samovar.MAG.mapping import validate_mapping_resume
    folder = tmp_path/'initial_annotations';folder.mkdir()
    (folder/'s.csv').write_text('seq,taxID_assembly_0,feat_assembly_0_MAP_MAG_ID\nr,562,M1\n')
    with pytest.raises(ValueError,match='Legacy MAP_MAG_ID'):
        validate_mapping_resume(tmp_path)


def test_nonempty_catalog_is_frozen_and_all_initial_samples_contribute(tmp_path, monkeypatch):
    from samovar.MAG.homology import stabilize_mags
    import numpy as np
    rng = np.random.default_rng(2)
    for sample, alphabet in (('s1','AACCGT'),('s2','ACGGTT')):
        source = tmp_path/sample;source.mkdir()
        sequence = ''.join(rng.choice(list(alphabet),3000))
        (source/'c.fa').write_text('>c\n'+sequence+'\n')
        stabilize_mags(source,tmp_path/'catalog.tsv',tmp_path/(sample+'_stable'))
    fastas = list((tmp_path/'fastas').glob('*.fa'))
    assert len(fastas) == 2
    def taxonomy(mags, dest, **kwargs):
        Path(dest).write_text('mag_id\ttaxid\tlineage\n' + ''.join(f'{p.stem}\t562\td__Bacteria\n' for p in fastas))
    monkeypatch.setattr('samovar.assembly_profiling.run_mag_taxonomy',taxonomy)
    finalize_mapping_catalog(tmp_path/'catalog.tsv', samples=['s2','s1'],taxonomy_tool='custom')
    validate_mapping_catalog(tmp_path/'catalog.tsv')
    (tmp_path/'fastas/extra.fa').write_text('>extra\nACGT\n')
    with pytest.raises(ValueError,match='FASTA set'):
        validate_mapping_catalog(tmp_path/'catalog.tsv')


def test_initial_and_regenerated_features_train_and_predict(tmp_path):
    from samovar.reprofiling import train_models, predict_taxid
    tax = tmp_path/'tax.tsv'
    tax.write_text('mag_id\ttaxid\tlineage\nM1\t562\t\nM2\t9606\t\n')
    frames = []
    for phase in ('initial','regenerated'):
        reads = fastq(tmp_path/(phase+'.fastq'), [f'r{i}' for i in range(40)])
        sam = tmp_path/(phase+'.sam')
        sam.write_text(''.join(alignment(f'r{i}', 'M1' if i%2 else 'M2', score=20) for i in range(40)))
        out = tmp_path/(phase+'.tsv')
        assign_alignments(reads,'',str(sam),str(tax),str(out))
        frame = pd.read_table(out).rename(columns={'taxID':'taxID_assembly_0','MAG_ID':'feat_assembly_0_MAG_ID'})
        frame['true'] = frame['taxID_assembly_0']
        frames.append(frame)
    initial, regenerated = frames
    model, _, _, features = train_models(regenerated, methods='random_forest',random_state=42)
    assert 'feat_assembly_0_MAG_ID' in features
    result = predict_taxid(initial.drop(columns='true'), model_path=model, feature_cols=features)
    assert list(result['taxid_SAMOVAR'].astype(int)) == list(initial['true'])


def test_empty_reference_maps_every_real_read_to_zero(tmp_path, monkeypatch):
    from samovar.assembly_profiling import run_assembly_annotator
    monkeypatch.setenv('SAMOVAR_MAG_CATALOG',str(tmp_path/'catalog.tsv'))
    monkeypatch.delenv('SAMOVAR_ASSEMBLY_EXTRA',raising=False)
    finalize_mapping_catalog(tmp_path/'catalog.tsv',samples=['s'])
    reads = fastq(tmp_path/'reads.fastq',['real'])
    monkeypatch.setattr('samovar.assembly_profiling.run_assembler',lambda *a,**k: pytest.fail('regenerated must not assemble'))
    out = run_assembly_annotator(reads,'',str(tmp_path/'out.tsv'),mag_mode='map-initial')
    assert out.read_text() == 'seq\ttaxID\tMAG_ID\nreal\t0\t\n'


def test_missing_marker_coordinates_is_corruption_not_biological_empty(tmp_path):
    (tmp_path/'gtdbtk.bac120.markers_summary.tsv').write_text(
        'name\tlist_unique_genes\tlist_multiple_genes\tlist_multiple_unique_genes\nM1\tPF1\t\t\n')
    with pytest.raises(ValueError,match='coordinates are missing'):
        extract_gtdb_markers(tmp_path,tmp_path/'markers.tsv')


def test_native_mapping_failure_is_not_zero_annotation(tmp_path, monkeypatch):
    import subprocess
    from samovar.MAG.mapping import write_mapping_table
    mags=tmp_path/'mags';mags.mkdir();(mags/'M1.fa').write_text('>c\nACGTACGTACGT\n')
    def fail(*args,**kwargs):
        raise subprocess.CalledProcessError(1,['minimap2'])
    monkeypatch.setattr('samovar.assembly_profiling.run_minimap2',fail)
    with pytest.raises(subprocess.CalledProcessError):
        write_mapping_table(fastq(tmp_path/'reads.fastq',['r']),'',str(mags),'',str(tmp_path/'out.tsv'))
    assert not (tmp_path/'out.tsv').exists()


def test_unclassified_mag_is_a_zero_result_not_invented_feature(tmp_path):
    table, out = resolve(tmp_path,alignment('unique','M1'),taxa=('0','0'))
    assert str(table.loc['unique','taxID']) == '0'
    assert table.loc['unique','MAG_ID'] == ''
    with sqlite3.connect(str(out)+'.candidates.sqlite') as db:
        assert db.execute('select count(*) from candidates').fetchone()[0] == 1


def test_empty_gtbd_summary_row_cannot_erase_assigned_taxon(tmp_path, monkeypatch):
    from samovar.assembly_profiling import run_gtdbtk
    root = tmp_path/'db'; (root/'metadata').mkdir(parents=True)
    monkeypatch.setenv('GTDBTK_DATA_PATH',str(root))
    monkeypatch.setattr('samovar.assembly_profiling.which_tool',lambda *a:'/configured/gtdbtk')
    monkeypatch.setattr('samovar.assembly_profiling._gtdb_taxid_from_lineage',lambda s: '562' if s else '0')
    def run(cmd, **kwargs):
        work = Path(cmd[cmd.index('--out_dir')+1]);work.mkdir()
        (work/'valid.summary.tsv').write_text('user_genome\tclassification\nM1\td__Bacteria;s__Escherichia coli\n')
        (work/'empty.summary.tsv').write_text('user_genome\tclassification\nM1\t\n')
    monkeypatch.setattr('samovar.assembly_profiling._run',run)
    out=run_gtdbtk(str(tmp_path/'mags'),str(tmp_path/'tax.tsv'))
    assert out.read_text().splitlines() == ['mag_id\ttaxid\tlineage','M1\t562\td__Bacteria;s__Escherichia coli']


def test_apply_maps_to_copied_reference_without_mutating_source_config(tmp_path):
    from samovar.apply import _retarget_init_config
    cfg = {'mag_catalog_root':'/source/MAG','mag_catalog_phase':'initial','run_config':[]}
    out = _retarget_init_config(cfg,tmp_path,2,512)
    assert out['mag_catalog_root'] == str(tmp_path/'MAG')
    assert out['mag_catalog_phase'] == 'regenerated'
    assert cfg['mag_catalog_root'] == '/source/MAG'


def test_ani_classified_reference_identifies_markers_without_reclassification(tmp_path,monkeypatch):
    from samovar.assembly_profiling import run_gtdbtk
    root=tmp_path/'db';(root/'metadata').mkdir(parents=True)
    monkeypatch.setenv('GTDBTK_DATA_PATH',str(root))
    monkeypatch.setattr('samovar.assembly_profiling.which_tool',lambda *a:'/configured/gtdbtk')
    monkeypatch.setattr('samovar.assembly_profiling._gtdb_taxid_from_lineage',lambda s:'562')
    calls=[]
    def run(cmd,**kwargs):
        calls.append(cmd)
        work=Path(cmd[cmd.index('--out_dir')+1]);work.mkdir(exist_ok=True)
        if cmd[1]=='classify_wf':
            (work/'valid.summary.tsv').write_text('user_genome\tclassification\nM1\td__Bacteria;s__Escherichia coli\n')
        else:
            assert cmd[1]=='identify'
            (work/'gtdbtk.bac120.markers_summary.tsv').write_text('name\tlist_unique_genes\nM1\t\n')
    monkeypatch.setattr('samovar.assembly_profiling._run',run)
    out=run_gtdbtk(str(tmp_path/'mags'),str(tmp_path/'tax.tsv'),require_markers=True)
    assert len(calls)==2
    assert '--keep_intermediates' in calls[0]
    assert 'M1\t562\td__Bacteria' in out.read_text()


def test_mapping_catalog_preserves_contig_boundaries(tmp_path):
    from samovar.MAG.homology import stabilize_mags
    from samovar.seqio import iter_seqio_fasta
    source=tmp_path/'mags';source.mkdir()
    (source/'bin.fa').write_text('>a\n'+'A'*1000+'\n>b\n'+'C'*1000+'\n')
    stabilize_mags(source,tmp_path/'catalog.tsv',tmp_path/'stable',preserve_contigs=True)
    fasta=next((tmp_path/'fastas').glob('*.fa'))
    records=list(iter_seqio_fasta(fasta))
    assert len(records)==2
    assert str(records[0].seq)=='A'*1000 and str(records[1].seq)=='C'*1000
    assert len({record.id for record in records})==2


def test_missing_query_is_not_successful_empty_output(tmp_path):
    with pytest.raises(FileNotFoundError):
        assign_alignments(str(tmp_path/'absent.fastq'),'','','',str(tmp_path/'out.tsv'))
    assert not (tmp_path/'out.tsv').exists()


def test_initial_cohort_cannot_silently_grow_after_freezing(tmp_path):
    finalize_mapping_catalog(tmp_path/'catalog.tsv',samples=['s1'])
    with pytest.raises(ValueError,match='cohort changed'):
        validate_mapping_catalog(tmp_path/'catalog.tsv',samples=['s1','s2'])


def test_native_mapping_honors_imported_candidate_flags(tmp_path,monkeypatch):
    from samovar.MAG.mapping import write_mapping_table
    config=tmp_path/'install.json'
    config.write_text(json.dumps({'tools':{'minimap2':{'exec':{'path':'/configured/minimap2'},'flags':'-N 7 --secondary=no'}}}))
    monkeypatch.setenv('SAMOVAR_CONFIG',str(config))
    mags=tmp_path/'mags';mags.mkdir();(mags/'M1.fa').write_text('>c\nACGTACGTACGT\n')
    tax=tmp_path/'tax.tsv';tax.write_text('mag_id\ttaxid\tlineage\nM1\t562\t\n')
    calls=[]
    def run(r1,r2,reference,bam,**kwargs):
        calls.append(kwargs['extra'])
        Path(bam).write_text(alignment('r','M1'))
    monkeypatch.setattr('samovar.assembly_profiling.run_minimap2',run)
    out=tmp_path/'out.tsv'
    write_mapping_table(fastq(tmp_path/'reads.fastq',['r']),'',str(mags),str(tax),str(out))
    assert calls==[['-N','7','--secondary=no','--secondary=yes']]
    with sqlite3.connect(str(out)+'.candidates.sqlite') as db:
        settings=json.loads(db.execute('select value from metadata where key="alignment_settings"').fetchone()[0])
        assert settings['flags']==calls[0]


def test_batch_mapping_tables_train_when_aggregate_csv_is_absent(tmp_path):
    from samovar.reprofilers import load_regenerated_table, run_reprofiler
    source=tmp_path/'regenerated';source.mkdir()
    frame=pd.DataFrame({'seq':[f'r{i}' for i in range(40)],
                        'taxID_assembly_0':[562 if i%2 else 9606 for i in range(40)],
                        'feat_assembly_0_MAG_ID':['M1' if i%2 else 'M2' for i in range(40)],
                        'true':[562 if i%2 else 9606 for i in range(40)]})
    frame.to_csv(source/'s.annotation.csv',index=False)
    loaded=load_regenerated_table(source/'combined_annotation_table.csv')
    assert len(loaded)==40
    out=tmp_path/'model'
    result=run_reprofiler('random_forest',regenerated_path=source/'combined_annotation_table.csv',
                         initial={'s':frame.drop(columns='true')},output_dir=out,
                         config={'cores':1,'feature_importance':'none'})
    assert result.model is not None
    assert (out/'trained_model.joblib').is_file()
