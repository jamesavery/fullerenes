"""Make `import fullerenes` resolve to the in-tree package.

Works whether or not `pip install -e .` has been run, as long as the extension
has been built into the source `fullerenes/` dir (the default cmake dev build).
"""
import os
import sys

_pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
