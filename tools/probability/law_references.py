#!/usr/bin/env python3
"""Valeurs de référence des lois de wsl::random, une famille de lois à la fois.

    uv run --no-project --with mpmath python tools/probability/law_references.py FAMILLE FICHIER

FAMILLE : un module de tools/probability/references/ (base, derived, rice…) ; CMake engendre
chaque famille dans le dossier de build et ne recalcule que celles dont le module change.
Non versionnées. La famille base vérifie aussi les constantes en double-double de
include/wsl/special_functions/detail/extended.hpp.
"""

import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "references"))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage : law_references.py FAMILLE FICHIER")
    family = sys.argv[1]
    sys.argv = [sys.argv[0], sys.argv[2]]
    importlib.import_module(family).main()
