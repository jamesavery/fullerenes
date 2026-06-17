"""Python bindings for the fullerene C++ library.

Topology, geometry, enumeration and canonical naming of fullerene isomers,
backed by libfullerenes.so. The binding surface is mostly generated from the
clang AST (see tools/gen_bindings.py); this package adds thin Python sugar and
sensible defaults on top of the compiled ``_fullerenes`` extension.
"""

import os as _os

from . import _fullerenes as _ext
from ._fullerenes import *  # noqa: F401,F403  (re-export the compiled surface)

__all__ = [n for n in dir(_ext) if not n.startswith("_")]
__version__ = getattr(_ext, "version", lambda: "0.0.0")()

# Database location. Honour the library's compiled-in FULLERENE_DATABASE_PATH
# (queryable via get_database_path()) when it exists; only when that path is
# absent fall back to the conventional local copy at ~/fullerene-data/database.
# Override at runtime with fullerenes.set_database_path(path) -- always a setter,
# never an environment variable. Only IsomerDB.read/get/make need the database;
# number_isomers/symmetries/buckygen are file-free.
_fallback_db = _os.path.expanduser("~/fullerene-data/database")
if not _os.path.isdir(_ext.get_database_path()) and _os.path.isdir(_fallback_db):
    _ext.set_database_path(_fallback_db)
