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
python3 -m pytest tests/ -v  # 57 tests
```

## Database

SQLite at `data/keloid.db`. 3 tables: `keloid_targets`, `drugs`, `drug_targets`. Join key: `gene_symbol` (HGNC standard).

## Scoring note

`OpenTargets_association` is excluded from pathway overlap count. Only curated pathways count toward `pathway_overlap_count` and the multi-target bonus.

## Commit style

Conventional commits, focused and atomic.
