"""Built-in grouped whole-row CAT resampler."""

from samovar.resample import resample_cat


def resample(initial, regenerated, groups, config):
    seed = 0 if config is None else int(config.get("seed", 0))
    columns = [] if config is None else config.get("columns", [])
    return resample_cat(initial, regenerated, groups, columns=columns, seed=seed)
