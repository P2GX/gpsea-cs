"""
GPSEA-CS is a repository with examples for GPSEA. This package is designed to structure the
summary output for each example notebook.
"""

__version__ = "0.0.1"


from .summarizer import GpseaSummarizer, SignificantResults
from .html_visualizer import HtmlVisualizer

__all__ = [
    "GpseaSummarizer",
    "HtmlVisualizer",
    "SignificantResults"
]
