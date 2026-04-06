#!/usr/bin/env bash
set -euo pipefail

echo "=== Step 1: Seeding keloid targets ==="
python3 -m scripts.01_seed_keloid_targets

echo "=== Step 2: Ingesting drugs from DGIdb ==="
python3 -m scripts.02_ingest_drugs

echo "=== Step 3: Scoring candidates ==="
python3 -m scripts.03_score_candidates

echo "=== Step 4: Cross-referencing clinical trials ==="
python3 -m scripts.05_clinical_trials

echo "=== Step 6: Structural docking (may take hours on first run) ==="
python3 -m scripts.06_structural_docking

echo "=== Step 5: Generating report ==="
python3 -m scripts.04_generate_report

echo "=== Done. Report at output/report.md ==="
