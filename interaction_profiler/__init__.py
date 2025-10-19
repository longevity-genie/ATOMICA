"""ATOMICA Interaction Profiler - Tools for analyzing protein-ligand interactions.

This module provides CLI tools for:
- Computing interaction scores (importance of blocks) using masked prediction scores
- Preparing molecular structures for visualization and analysis
"""

from .interact_score import app

__all__ = ["app"]
