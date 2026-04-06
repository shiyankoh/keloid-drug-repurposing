# CLAUDE.md

Instructions for AI assistants working on this codebase.

## Module pattern

Each pipeline step has two files:
- **Core module** (`scripts/foo.py`) — importable functions, used by tests
- **Entry point** (`scripts/0N_foo.py`) — thin wrapper calling `main()`

Tests import from the core module. Pipeline runner uses `python3 -m scripts.0N_foo`.

## Running

```bash
source venv/bin/activate       # activate venv first
./run_pipeline.sh              # full pipeline (stages 01-06)
python3 -m pytest tests/ -v   # 118 tests
```

Vina binary is at `~/bin/vina`. Make sure `$HOME/bin` is in PATH.

## Database

SQLite at `data/keloid.db`. 6 tables:

| Table | Purpose |
|-------|---------|
| `keloid_targets` | 22 curated + ~500 OpenTargets gene targets |
| `drugs` | ~1,230 FDA-approved drugs from DGIdb |
| `drug_targets` | Drug-gene interactions. Join key: `gene_symbol` |
| `protein_structures` | PDB/AlphaFold 3D structures for dockable targets |
| `drug_structures` | PubChem SDF + PDBQT files for top drugs |
| `docking_results` | AutoDock Vina binding energies per drug-protein pair |

## Pipeline stages

1. `seed_keloid_targets` — curated CSV + OpenTargets API
2. `ingest_drugs` — DGIdb API for FDA-approved drugs
3. `score_candidates` — composite score (pathway overlap + evidence + multi-target)
4. `generate_report` — Markdown report with severity tiers
5. `clinical_trials` — ClinicalTrials.gov cross-reference
6. `06_structural_docking` — fetch structures, prep proteins, run Vina docking

## Structural docking (Stage 06)

Modules:
- `scripts/protein_structures.py` — fetch PDB/AlphaFold protein structures
- `scripts/drug_structures.py` — fetch PubChem drug structures, convert to PDBQT
- `scripts/protein_prep.py` — clean PDB files, extract binding sites from HETATM
- `scripts/docking.py` — run AutoDock Vina, parse output, store results

Dependencies: `brew install open-babel`, Vina binary at `~/bin/vina`.

Docking is CPU-intensive (~2-10 min per pair). Full batch of ~20 pairs takes ~30 min. Structure files are cached in `data/structures/` (gitignored). API responses cached in `data/raw/`.

Only 12 of 22 curated targets are dockable (kinases/enzymes). Structural proteins (COL1A1/A2, COL3A1) and secreted ligands (TGFB1/B2, WNT3A, VEGFA, IL6) are skipped. Biologics (rituximab, aflibercept, etc.) can't be fetched from PubChem — they're large proteins, not small molecules.

## Scoring

Two scoring modes:
- **Base score:** `pathway_overlap × 0.40 + evidence × 0.40 + multi_target × 0.20`
- **With structural data:** `overlap × 0.30 + evidence × 0.30 + structural × 0.25 + multi_target × 0.15`

`OpenTargets_association` is excluded from pathway overlap count. Only curated pathways count toward `pathway_overlap_count` and the multi-target bonus.

Binding energy → structural score: -10 kcal/mol → 1.0, -3 kcal/mol → 0.0 (linear scale, clamped).

## Commit style

Conventional commits, focused and atomic.
