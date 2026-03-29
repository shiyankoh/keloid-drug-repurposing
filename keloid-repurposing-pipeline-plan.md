# Keloid Drug Repurposing Data Pipeline

A personal research pipeline to identify promising repurposed drug candidates for keloids, modeled on Fajgenbaum's computational pharmacophenomics approach.

## Core Thesis

Keloids are driven by well-characterized molecular pathways (TGF-β, mTOR, PI3K/AKT, IL-6, collagen synthesis) that overlap with pathways targeted by drugs approved for other diseases. A data pipeline can systematically surface these overlaps and rank candidates by evidence strength.

---

## v1 Scope (Reduced)

v1 focuses on the core overlap query: which FDA-approved drugs hit keloid-relevant molecular targets? Everything else (literature search, clinical trials, safety profiles) is deferred to follow-up phases. See TODOS.md for deferred work.

### Architecture

```
                        ┌─────────────────┐
                        │  Python Scripts  │
                        │  (data ingest)   │
                        └────────┬────────┘
                                 │
              ┌──────────────────┼──────────────────┐
              ▼                  ▼                   ▼
     ┌────────────────┐ ┌───────────────┐  ┌────────────────┐
     │  OpenTargets    │ │   DGIdb       │  │  Manual seed   │
     │  GraphQL API    │ │   REST API    │  │  (curated CSV) │
     │  (targets +     │ │  (drug-gene   │  │  (pathways)    │
     │   drug assoc.)  │ │   links)      │  │                │
     └────────┬────────┘ └───────┬───────┘  └────────┬───────┘
              │                  │                    │
              ▼                  ▼                    ▼
    ┌───────────────────────────────────────────────────────────┐
    │                   Raw JSON cache                          │
    │                   (data/raw/)                             │
    └──────────────────────┬───────────────────────────────────┘
                           │
                           ▼
        ┌──────────────────────────────────────────────────┐
        │              SQLite (local)                      │
        │                                                  │
        │  keloid_targets  ←──→  drug_targets  ←──→ drugs  │
        │                                                  │
        └──────────────────────┬───────────────────────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │  Scoring query (SQL) │
                    │  + Python post-      │
                    │    processing        │
                    └──────────┬──────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │  Markdown report    │
                    │  (ranked candidates)│
                    └─────────────────────┘
```

### Project Structure

```
keloid/
├── scripts/
│   ├── 01_seed_keloid_targets.py   # Curated seed + OpenTargets pull
│   ├── 02_ingest_drugs.py          # DGIdb + OpenTargets drug data
│   ├── 03_score_candidates.py      # Core overlap query + ranking
│   └── 04_generate_report.py       # Markdown output
├── tests/
│   ├── test_normalization.py       # Gene symbol normalization
│   ├── test_dedup.py               # Deduplication logic
│   ├── test_scoring.py             # Score calculation + edge cases
│   └── fixtures/                   # Sample API responses for mocking
├── data/
│   ├── raw/                        # Cached API responses (JSON)
│   ├── keloid.db                   # SQLite database (generated)
│   └── seed_targets.csv            # Manual keloid pathway seeds
├── output/
│   └── report.md                   # Generated candidate report
├── run_pipeline.sh                 # Runs all 4 scripts in order
├── requirements.txt
├── CLAUDE.md
└── TODOS.md
```

---

## Data Model (SQLite)

```
keloid_targets
├── id (integer, PK, autoincrement)
├── target_name (text)
├── gene_symbol (text, UNIQUE)
├── ensembl_id (text)               # For OpenTargets cross-reference
├── pathway (text)
├── evidence_type (text)            # genetic | expression | functional
├── evidence_strength (real)        # 0.0–1.0 (OpenTargets score or 0.7 default)
├── source (text)                   # pmid or database name
└── created_at (text)               # ISO 8601

drugs
├── id (integer, PK, autoincrement)
├── drug_name (text)
├── generic_name (text, UNIQUE)
├── approval_status (text)
├── original_indication (text)
├── mechanism_of_action (text)
└── created_at (text)

drug_targets
├── id (integer, PK, autoincrement)
├── drug_id (integer → drugs)
├── gene_symbol (text)              # Joins to keloid_targets.gene_symbol
├── action_type (text)              # inhibitor | agonist | modulator
├── source (text)                   # DGIdb, DrugBank, etc.
└── created_at (text)
```

---

## Phase 1: Build the Knowledge Base

### 1A. Keloid Pathway Map

Pull and structure data on every known keloid-associated pathway and molecular target.

**Sources:** OpenTargets GraphQL API + manual seed list (CSV)

**Key pathways to seed manually:**
- TGF-β/Smad
- PI3K/AKT/mTOR
- Wnt/β-catenin
- JAK/STAT
- MAPK/ERK
- EGFR
- IL-6/TNF-α signaling
- Collagen biosynthesis (PRS/proline)
- CYP24A1/vitamin D metabolism
- VEGF (angiogenesis)

**Evidence scoring:**
- OpenTargets targets: use their association score directly (0.0–1.0)
- Manually seeded targets: default to 0.7

**Gene ID normalization:**
- OpenTargets returns Ensembl gene IDs (e.g., `ENSG00000198793`)
- Store both `ensembl_id` and `gene_symbol` (HGNC)
- Normalization step during ingest maps Ensembl → HGNC using OpenTargets' own mapping
- DGIdb returns HGNC gene symbols natively — no mapping needed

### 1B. Drug-Target Database

Map FDA-approved drugs to their known molecular targets.

**Source:** DGIdb REST API (aggregates DrugBank, ChEMBL, and others)

**Process:**
1. Query DGIdb for all interactions with genes in `keloid_targets`
2. Filter to FDA-approved drugs
3. Split into `drugs` + `drug_targets` records
4. Deduplicate drugs by `generic_name`

**API caching:** All raw API responses saved to `data/raw/` as JSON before processing. Scripts can re-run from cache without hitting the network.

---

## Phase 2: Cross-Reference & Scoring

### 2A. Pathway Overlap Scoring

The core query: which approved drugs hit keloid-relevant targets?

```sql
SELECT d.drug_name, d.original_indication, dt.gene_symbol, kt.pathway,
       kt.evidence_strength
FROM drugs d
JOIN drug_targets dt ON d.id = dt.drug_id
JOIN keloid_targets kt ON dt.gene_symbol = kt.gene_symbol
WHERE d.approval_status = 'approved'
ORDER BY kt.evidence_strength DESC;
```

### 2B. Composite Scoring (v1)

```
v1_score = (
    pathway_overlap_count × 0.40 +
    avg_evidence_strength  × 0.40 +
    multi_target_bonus     × 0.20
)
```

- **pathway_overlap_count:** Number of distinct keloid pathways the drug touches
- **avg_evidence_strength:** Average evidence_strength across matched targets
- **multi_target_bonus:** Drugs hitting 3+ distinct pathways get a bonus (Fajgenbaum's insight — diseases are multi-pathway, so multi-target drugs may be more effective)

### 2C. Output: Ranked Candidate Report

For each top candidate, generate a Markdown brief containing:
- Drug name, current approved indication, mechanism of action
- Which keloid pathways it targets and supporting evidence
- Number of pathway overlaps and composite score
- Suggested next step (discuss with dermatologist, literature deep-dive, etc.)

---

## Technical Stack

| Component | Tool | Notes |
|-----------|------|-------|
| Database | SQLite (local) | Zero-config, portable, sufficient for 3 small tables |
| Data ingestion | Python scripts | `requests` for APIs |
| Orchestration | Standalone scripts + `run_pipeline.sh` | Run individually or all at once |
| Testing | pytest | Unit tests on transforms + scoring logic |
| Output | Markdown report | Human-readable, portable |

## API Endpoints

| Source | API | Auth | Rate Limits |
|--------|-----|------|-------------|
| OpenTargets | api.platform.opentargets.org (GraphQL) | None | Open |
| DGIdb | dgidb.org/api | None | Open |

---

## Key Design Decisions

| # | Decision | Rationale |
|---|----------|-----------|
| 1 | gene_symbol as join key + ensembl_id column | DGIdb returns HGNC symbols natively. Normalize OpenTargets Ensembl IDs during ingest. Explicit > clever. |
| 2 | OpenTargets scores + 0.7 default for manual seeds | Consistent, evidence-based. Good enough for a hypothesis generator. |
| 3 | Local SQLite (not Supabase) | Zero-config, no credentials, no network dependency. Migrate to Supabase only if/when pgvector or a dashboard is needed. |
| 4 | Standalone scripts + run_pipeline.sh | Maximum visibility and debuggability. Inspect output between steps. |
| 5 | Cache raw API responses to JSON | Insurance against flaky academic APIs. Re-runnable without network. |
| 6 | pytest unit tests on transforms + scoring | Test the logic, not the plumbing. ~30-40 tests covering normalization, dedup, scoring, edge cases. |

---

## Known Critical Gaps (must address during implementation)

1. **OpenTargets returns 0 keloid targets** — must log a warning and fail gracefully, not produce an empty report silently
2. **Ensembl ID with no HGNC mapping** — must log a warning with the unmapped ID, not silently drop the target

---

## Quick Wins (do before or alongside the pipeline)

1. **Search Every Cure's MeDIC database** (open-source) for keloid or fibrosis entries
2. **Query OpenTargets web UI** for "keloid" — see which drugs already have association scores
3. **PubMed search:** `"drug repurposing" AND keloid` — see what's been published
4. **Ask your dermatologist** about their experience with off-label treatments (verapamil, 5-FU, bleomycin, sirolimus)
5. **Check if Every Cure accepts disease nominations** — keloids affect millions with no gold-standard therapy

---

## Important Caveats

- This pipeline generates **hypotheses**, not medical advice. Any candidate must be discussed with a dermatologist before use.
- Off-label prescribing is legal and common, but requires a willing physician who understands the evidence.
- Keloid treatment response varies significantly by genetics, location, and scar maturity — a drug that works mechanistically may still not work for a specific case.
- Fajgenbaum's team has PhDs in bioinformatics and access to medical records data you won't have — this is a simplified version focused on publicly available data.
