# CLAUDE.md

Instructions for AI assistants working on this codebase.

## Module pattern

Each pipeline step has two files:
- **Core module** (`scripts/foo.py`) — importable functions, used by tests
- **Entry point** (`scripts/0N_foo.py`) — thin wrapper calling `main()`

Tests import from the core module. Pipeline runner uses `python3 -m scripts.0N_foo`.

## Running

```bash
./run_pipeline.sh          # full pipeline
source venv/bin/activate      # activate venv first
python3 -m pytest tests/ -v  # 118 tests
```

## Database

SQLite at `data/keloid.db`. 6 tables: `keloid_targets`, `drugs`, `drug_targets`, `protein_structures`, `drug_structures`, `docking_results`. Join key: `gene_symbol` (HGNC standard).

## Structural Docking (Stage 06)

New modules for protein-drug molecular docking:
- `scripts/protein_structures.py` — fetch PDB/AlphaFold protein structures
- `scripts/drug_structures.py` — fetch PubChem drug structures, convert to PDBQT
- `scripts/protein_prep.py` — clean PDB files, find binding sites
- `scripts/docking.py` — run AutoDock Vina, parse results

Requires: `brew install open-babel`, Vina binary at `~/bin/vina`, Python venv with `requests` + `pytest`.

Docking is CPU-intensive (~2-10 min per pair). Full batch of ~60 pairs takes hours.

## Scoring note

`OpenTargets_association` is excluded from pathway overlap count. Only curated pathways count toward `pathway_overlap_count` and the multi-target bonus.

## Commit style

Conventional commits, focused and atomic.
