"""Shared cohort clustering is reusable across annotation phases."""
import pandas as pd
import pytest

from samovar.kmer_cluster import predict, train, validate_model
from samovar.reprofiling import preprocess_data


def fastq(path, entries):
    path.write_text(''.join(f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n' for name, seq in entries))
    return str(path)


def test_shared_model_predict_and_categorical_schema(tmp_path):
    a = fastq(tmp_path / 'a.fastq', [('a', 'AAAAAA'), ('b', 'CCCCCC')])
    b = fastq(tmp_path / 'b.fastq', [('c', 'GGGGGG'), ('d', 'TTTTTT')])
    model = tmp_path / 'model.joblib'
    trained = train([(a, ''), (b, '')], model, k=2, min_clusters=2, max_clusters=2, seed=7)
    assert trained['training_reads'] == 4
    predict(a, '', model, tmp_path / 'initial.tsv')
    predict(a, '', model, tmp_path / 'regenerated.tsv')
    initial = pd.read_csv(tmp_path / 'initial.tsv', sep='\t')
    regenerated = pd.read_csv(tmp_path / 'regenerated.tsv', sep='\t')
    pd.testing.assert_frame_equal(initial, regenerated)
    # Cluster identifiers preserve their numeric identity even if one cluster
    # is absent from the reprofiler's training data.
    col = 'feat_cluster_0_kmer_cluster'
    raw = pd.DataFrame({col: [0, 0], 'true': [1, 1]})
    raw.attrs["samovar_cluster_levels"] = {col: [0, 1]}
    fit = preprocess_data(raw)
    schema = fit.attrs['samovar_preprocessing']
    got = preprocess_data(pd.DataFrame({col: [1], 'true': [1]}), preprocessing=schema)
    assert got[f"{col}__cluster_1"].iloc[0] == 1
    assert got[f"{col}__cluster_0"].iloc[0] == 0
    validate_model(model, '-k 2 --min-clusters 2 --max-clusters 2 --seed 7', [(a, ''), (b, '')])
    with pytest.raises(ValueError, match='settings'):
        validate_model(model, '-k 3')


def test_balanced_training_paired_and_empty(tmp_path):
    r1 = fastq(tmp_path / 'r1.fastq', [('a/1', 'AAAA'), ('b/1', 'CCCC')])
    r2 = fastq(tmp_path / 'r2.fastq', [('a/2', 'TTTT'), ('b/2', 'GGGG')])
    other = fastq(tmp_path / 'other.fastq', [('c', 'ATAT')])
    model = tmp_path / 'model.joblib'
    train([(r1, r2), (other, '')], model, k=2, max_clusters=2, sampling='balanced', max_reads=2)
    predict(r1, r2, model, tmp_path / 'out.tsv')
    assert len(pd.read_csv(tmp_path / 'out.tsv', sep='\t')) == 2
    invalid = fastq(tmp_path / 'invalid.fastq', [('x', 'NNNN')])
    with pytest.raises(ValueError, match='no reads'):
        train([(invalid, '')], tmp_path / 'invalid.joblib', k=2)
    predict(invalid, '', model, tmp_path / 'invalid.tsv')
    assert pd.read_csv(tmp_path / 'invalid.tsv', sep='\t').iloc[0, 1] != pd.read_csv(tmp_path / 'invalid.tsv', sep='\t').iloc[0, 1]


def test_missing_model_blocks_resume(tmp_path):
    import yaml
    from samovar.exec_control import check_startpoint
    config = tmp_path / '.log' / 'configs'
    config.mkdir(parents=True)
    (config / 'config_init.yaml').write_text(yaml.safe_dump({'run_config': [
        {'run_name': 'clusters', 'type': 'kmer_cluster', 'db_path': '.'}]}))
    markers = tmp_path / '.log' / 'checkpoints'
    markers.mkdir()
    (markers / 'annotate_initial.done').touch()
    assert any('saved kmer-cluster model' in gap for gap in check_startpoint(tmp_path))


def test_reverse_complements_share_cluster_and_mode_is_saved(tmp_path):
    reads = fastq(tmp_path / "strand.fastq", [("a", "AAAAAA"), ("t", "TTTTTT"), ("c", "CCCCCC")])
    model = tmp_path / "canonical.joblib"
    trained = train([(reads, "")], model, k=2, min_clusters=2, max_clusters=2)
    assert len(trained["kmers"]) == 10
    predict(reads, "", model, tmp_path / "features.tsv")
    frame = pd.read_csv(tmp_path / "features.tsv", sep="\t").set_index("seq")
    assert frame.loc["a", "feat_kmer_cluster"] == frame.loc["t", "feat_kmer_cluster"]
    strand = train([(reads, "")], tmp_path / "strand.joblib", k=2, max_clusters=2,
                   reverse_complemented=False)
    assert len(strand["kmers"]) == 16
    assert not strand["reverse_complemented"]


def test_full_cluster_vocabulary_is_loaded_from_run(tmp_path):
    import yaml
    from samovar.kmer_cluster import cluster_levels_for_run
    reads = fastq(tmp_path / 'train.fastq', [('a', 'AAAAAA'), ('c', 'CCCCCC')])
    config = tmp_path / '.log' / 'configs'
    config.mkdir(parents=True)
    (config / 'config_init.yaml').write_text(yaml.safe_dump({'run_config': [
        {'run_name': 'dummy', 'type': 'constant_taxid'},
        {'run_name': 'cluster', 'type': 'kmer_cluster', 'db_path': '.'}]}))
    path = tmp_path / 'kmer_cluster_models' / 'cluster' / 'model.joblib'
    train([(reads, '')], path, k=2, min_clusters=2, max_clusters=2)
    col = 'feat_kmer_cluster_1_kmer_cluster'
    assert cluster_levels_for_run(tmp_path, [col]) == {col: [0, 1]}
