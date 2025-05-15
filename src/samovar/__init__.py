"""
Python wrapper for the samovaR package.
"""

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

from .core import *  # noqa 