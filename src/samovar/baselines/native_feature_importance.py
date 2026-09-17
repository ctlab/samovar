"""Built-in feature-importance scorer (native attributes, then permutation)."""

from samovar.feature_importance import score_feature_importance as _score


def score_feature_importance(model, annotation, initial_abundance, regenerated_abundance, config):
    return _score(model, annotation, initial_abundance, regenerated_abundance, config)
