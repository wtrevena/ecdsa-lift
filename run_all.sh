#!/usr/bin/env bash
# Run every experiment in order. Each script writes its own JSON into
# results/ and prints a short summary.
set -e
cd "$(dirname "$0")"
PY="${PY:-.venv/bin/python}"
if [ ! -x "$PY" ]; then PY=python3; fi

echo "==> phase 1 (Teichmüller δ injectivity)"
"$PY" experiments/phase1_injectivity.py
echo
echo "==> phase 2 (δ structure probe)"
"$PY" experiments/phase2_structure.py
echo
echo "==> phase 4 (cohomological splitting)"
"$PY" experiments/phase4_homomorphic_section.py
