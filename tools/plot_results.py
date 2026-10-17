#!/usr/bin/env python3
"""Kept for convenience; same as the `ekpn-viz` command (pip install -e .)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from ekpn_viz.cli import main  # noqa: E402

if __name__ == "__main__":
    main()
