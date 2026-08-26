"""
Automated benchmarcing with SamovaR
"""

from .version import get_version

__version__ = get_version()

from .fasta_processor import * 
from .build_database import *
from .genome_fetcher import *