"""
Automated benchmarcing with SamovaR
"""

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

#from .core import *  
from .fasta_processor import * 
from .build_database import * 