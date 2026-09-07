#!/usr/bin/env python3
"""Compatibility wrapper for the cleaned random-forest workflow."""

from pathlib import Path
import runpy

script = Path(__file__).with_name("random_forest.py")
runpy.run_path(str(script), run_name="__main__")

