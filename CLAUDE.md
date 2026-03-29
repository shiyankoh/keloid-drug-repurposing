# Keloid Drug Repurposing Pipeline

## What This Is

Personal research pipeline to find FDA-approved drugs that target keloid-relevant molecular pathways. Batch data pipeline — no web app. Outputs a ranked Markdown report with severity tiers and clinical trial cross-references.

## Project Structure

```
keloid/
├── scripts/
│   ├── db.py                       # Shared DB layer (SQLite init + connection)
│   ├── seed_keloid_targets.py      # Core module: CSV + OpenTargets → keloid_targets
│   ├── 01_seed_keloid_targets.py   # Entry point wrapper
│   ├── ingest_drugs.py             # Core module: DGIdb → drugs + drug_targets
│   ├── 02_ingest_drugs.py          # Entry point wrapper
│   ├── score_candidates.py         # Core module: overlap join + scoring
│   ├── 03_score_candidates.py      # Entry point wrapper
│   ├── clinical_trials.py          # Core module: ClinicalTrials.gov cross-ref
│   ├── 05_clinical_trials.py       # Entry point wrapper
│   ├── generate_report.py          # Core module: Markdown report generation
│   └── 04_generate_report.py       # Entry point wrapper
├── tests/
│   ├── conftest.py                 # Shared fixture: in-memory SQLite
│   ├── test_db.py
│   ├── test_seed_targets.py
│   ├── test_ingest_drugs.py
│   ├── test_scoring.py
│   ├── test_clinical_trials.py
│   └── fixtures/                   # Sample API responses for mocking
├── data/
│   ├── raw/                        # Cached API responses (JSON) — gitignored
│   ├── keloid.db                   # SQLite database — gitignored
│   ├── seed_targets.csv            # 22 manually curated keloid pathway targets
│   └── severity_tiers.csv          # Safety classification for top 50 drugs
├── output/                         # Generated reports — gitignored
├── run_pipeline.sh                 # Runs all 5 steps in sequence
├── requirements.txt
├── CLAUDE.md
└── TODOS.md
```

## Module Pattern

Each pipeline step has two files:
- **Core module** (`scripts/foo.py`) — importable functions, used by tests
- **Entry point** (`scripts/0N_foo.py`) — thin wrapper calling `main()`

Tests import from the core module. Pipeline runner uses `python3 -m scripts.0N_foo`.

## How to Run

```bash
# Full pipeline (5 steps)
./run_pipeline.sh

# Individual steps (must use -m flag for module imports)
python3 -m scripts.01_seed_keloid_targets
python3 -m scripts.02_ingest_drugs
python3 -m scripts.03_score_candidates
python3 -m scripts.05_clinical_trials
python3 -m scripts.04_generate_report

# Tests (57 tests)
python3 -m pytest tests/ -v
```

## Tech Stack

- **Python 3** — scripts + tests
- **SQLite** — local database, no credentials needed
- **pytest** — 57 unit tests
- **APIs:** OpenTargets (GraphQL, no auth), DGIdb (GraphQL, no auth), ClinicalTrials.gov (REST, no auth)

## Database Schema

3 tables in SQLite (`data/keloid.db`):

- `keloid_targets` — keloid-associated molecular targets with pathways and evidence scores
- `drugs` — FDA-approved drugs with indications and mechanisms
- `drug_targets` — junction table linking drugs to gene targets

Join key: `gene_symbol` (HGNC standard).

## Scoring

```
v1_score = (
    normalized_overlap   × 0.40    # curated pathway count / max
    avg_evidence_strength × 0.40    # mean across matched targets
    multi_target_bonus    × 0.20    # 1.0 if 3+ curated pathways
)
```

**Important:** `OpenTargets_association` is excluded from pathway overlap count. Only curated pathways (TGF-beta, VEGF, MAPK/ERK, etc.) count toward pathway_overlap_count and the multi-target bonus. OpenTargets genes still contribute to avg_evidence_strength.

## Severity Tiers

`data/severity_tiers.csv` classifies the top 50 drugs:
- **green** — feasible for keloid use (OTC, topical, well-known safety)
- **yellow** — moderate risk, needs evaluation (immunomodulators, repurposed drugs)
- **red** — oncology-grade toxicity, not appropriate for a benign condition
- **unclassified** — not yet reviewed

Report leads with "Actionable Candidates" — green/yellow drugs with 2+ curated pathways.

## Clinical Trials Cross-Reference

Step 5 fetches all keloid + hypertrophic scar trials from ClinicalTrials.gov and matches them to candidates using drug name aliases (e.g., sirolimus ↔ rapamycin). Results appear in the report with clickable NCT links.

## Key Design Decisions

1. **SQLite, not Supabase** — zero-config, no network dependency
2. **gene_symbol as join key** — DGIdb returns HGNC natively. OpenTargets normalized during ingest.
3. **Evidence scores:** OpenTargets association scores used directly. Manual seeds default to 0.7.
4. **Raw API responses cached** to `data/raw/` as JSON. Delete cache to re-fetch.
5. **OpenTargets_association excluded from pathway count** — prevents low-confidence noise from inflating scores

## Known Limitations

- Pipeline can't distinguish direct inhibition from indirect/downstream effects. A drug that "touches" a pathway may not therapeutically modulate it.
- DGIdb interaction types are often "unknown" — the drug interacts with the gene but the nature (inhibition vs. activation) isn't specified.
- Scoring rewards pathway breadth. A highly targeted drug (e.g., verapamil on calcium channels) can be clinically effective through a single mechanism but will score low.
- Clinical trial results are not analyzed — only existence of trials is captured.

## What's NOT Built Yet

See `TODOS.md`:
- Literature corpus + pgvector semantic search (needs Supabase migration)
- FDA FAERS safety profiles (partially addressed by severity tiers)
- Every Cure MeDIC database check (manual task)

## Commit Style

Conventional commits, focused and atomic.
