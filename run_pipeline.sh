#!/usr/bin/env bash
set -euo pipefail

echo "=== Step 1: Seeding keloid targets ==="
python3 scripts/01_seed_keloid_targets.py

echo "=== Step 2: Ingesting drugs from DGIdb ==="
python3 scripts/02_ingest_drugs.py

echo "=== Step 3: Scoring candidates ==="
python3 scripts/03_score_candidates.py

echo "=== Step 4: Generating report ==="
python3 scripts/04_generate_report.py

echo "=== Done. Report at output/report.md ==="
