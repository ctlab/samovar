"""
Automated benchmarcing with SamovaR
"""

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "0.10.6"

from .fasta_processor import * 
from .build_database import *
from .genome_fetcher import *