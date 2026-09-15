"""Built-in annotation converter: dump abundance tables; load them back."""

from samovar.abundance import load_table_input
from samovar.parse_annotators import Annotation


def dump(annotation: Annotation, dest, config):
    _ = config
    return annotation.write(dest, "abundance")


def load(path, config):
    _ = config
    return load_table_input(path)
