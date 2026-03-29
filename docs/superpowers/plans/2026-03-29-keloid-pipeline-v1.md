# Keloid Drug Repurposing Pipeline v1 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a batch data pipeline that surfaces FDA-approved drugs targeting keloid-relevant molecular pathways, ranked by evidence strength, output as a Markdown report.

**Architecture:** 4 standalone Python scripts run in sequence. Script 01 seeds keloid targets from OpenTargets + a curated CSV. Script 02 ingests drug-gene interactions from DGIdb. Script 03 joins and scores. Script 04 generates the Markdown report. All data lives in a local SQLite database. Raw API responses cached to JSON.

**Tech Stack:** Python 3, SQLite, pytest, requests. APIs: OpenTargets GraphQL (no auth), DGIdb GraphQL (no auth).

---

## File Structure

```
keloid/
├── scripts/
│   ├── db.py                       # Shared: DB init, connection, schema
│   ├── 01_seed_keloid_targets.py   # OpenTargets fetch + CSV merge → keloid_targets
│   ├── 02_ingest_drugs.py          # DGIdb fetch → drugs + drug_targets
│   ├── 03_score_candidates.py      # Overlap join + scoring → stdout summary
│   └── 04_generate_report.py       # Full Markdown report → output/report.md
├── tests/
│   ├── conftest.py                 # Shared fixtures: in-memory SQLite, sample data
│   ├── test_db.py                  # Schema creation tests
│   ├── test_seed_targets.py        # Normalization, merge, dedup, edge cases
│   ├── test_ingest_drugs.py        # DGIdb parsing, filtering, dedup
│   ├── test_scoring.py             # Score calc, multi-target bonus, ties, zero overlap
│   └── fixtures/
│       ├── opentargets_keloid.json # Sample OpenTargets response (3-5 targets)
│       └── dgidb_interactions.json # Sample DGIdb response (5-10 interactions)
├── data/
│   ├── raw/                        # Cached API responses (created at runtime)
│   ├── seed_targets.csv            # Manually curated keloid pathway seeds
│   └── keloid.db                   # SQLite DB (created at runtime, gitignored)
├── output/                         # Generated reports (gitignored)
├── run_pipeline.sh
├── requirements.txt
├── .gitignore
├── CLAUDE.md
└── TODOS.md
```

**Key boundary:** `scripts/db.py` is the only file that touches SQLite directly. All other scripts import from it. This keeps the DB layer in one place and makes testing easy (swap in an in-memory DB).

---

## Task 0: Project Setup

**Files:**
- Create: `requirements.txt`
- Create: `.gitignore`
- Create: `data/seed_targets.csv`
- Create: `run_pipeline.sh`

- [ ] **Step 1: Initialize git repo**

```bash
cd ~/Coding\ Projects/Keloid
git init
```

- [ ] **Step 2: Create requirements.txt**

```
requests>=2.31.0
pytest>=7.4.0
```

- [ ] **Step 3: Create .gitignore**

```
data/keloid.db
data/raw/
output/
__pycache__/
*.pyc
.pytest_cache/
venv/
.env
```

- [ ] **Step 4: Create the seed targets CSV**

This is the manually curated list of keloid-associated genes and pathways. Each row is one gene known to be involved in keloid pathogenesis. `evidence_strength` is 0.7 (our default for manually seeded targets — see plan decision #2).

```csv
gene_symbol,target_name,pathway,evidence_type,evidence_strength,source
TGFB1,transforming growth factor beta 1,TGF-β/Smad,functional,0.7,manual_curation
TGFB2,transforming growth factor beta 2,TGF-β/Smad,functional,0.7,manual_curation
SMAD3,SMAD family member 3,TGF-β/Smad,genetic,0.7,manual_curation
SMAD2,SMAD family member 2,TGF-β/Smad,functional,0.7,manual_curation
MTOR,mechanistic target of rapamycin kinase,PI3K/AKT/mTOR,functional,0.7,manual_curation
PIK3CA,phosphatidylinositol-4_5-bisphosphate 3-kinase catalytic subunit alpha,PI3K/AKT/mTOR,functional,0.7,manual_curation
AKT1,AKT serine/threonine kinase 1,PI3K/AKT/mTOR,functional,0.7,manual_curation
CTNNB1,catenin beta 1,Wnt/β-catenin,functional,0.7,manual_curation
WNT3A,Wnt family member 3A,Wnt/β-catenin,expression,0.7,manual_curation
JAK2,Janus kinase 2,JAK/STAT,functional,0.7,manual_curation
STAT3,signal transducer and activator of transcription 3,JAK/STAT,expression,0.7,manual_curation
MAPK1,mitogen-activated protein kinase 1,MAPK/ERK,functional,0.7,manual_curation
MAPK3,mitogen-activated protein kinase 3,MAPK/ERK,functional,0.7,manual_curation
EGFR,epidermal growth factor receptor,EGFR,expression,0.7,manual_curation
IL6,interleukin 6,IL-6/TNF-α,expression,0.7,manual_curation
TNF,tumor necrosis factor,IL-6/TNF-α,functional,0.7,manual_curation
COL1A1,collagen type I alpha 1 chain,Collagen biosynthesis,expression,0.7,manual_curation
COL1A2,collagen type I alpha 2 chain,Collagen biosynthesis,expression,0.7,manual_curation
COL3A1,collagen type III alpha 1 chain,Collagen biosynthesis,expression,0.7,manual_curation
CYP24A1,cytochrome P450 family 24 subfamily A member 1,Vitamin D metabolism,genetic,0.7,manual_curation
VEGFA,vascular endothelial growth factor A,VEGF,expression,0.7,manual_curation
KDR,kinase insert domain receptor,VEGF,functional,0.7,manual_curation
```

- [ ] **Step 5: Create run_pipeline.sh**

```bash
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
```

- [ ] **Step 6: Create directory structure**

```bash
mkdir -p scripts tests/fixtures data/raw output
chmod +x run_pipeline.sh
```

- [ ] **Step 7: Install dependencies and commit**

```bash
python3 -m pip install -r requirements.txt
git add requirements.txt .gitignore data/seed_targets.csv run_pipeline.sh CLAUDE.md TODOS.md keloid-repurposing-pipeline-plan.md docs/
git commit -m "chore: project setup with seed data and pipeline runner"
```

---

## Task 1: Database Layer (`scripts/db.py`)

**Files:**
- Create: `scripts/db.py`
- Create: `tests/conftest.py`
- Create: `tests/test_db.py`

- [ ] **Step 1: Write the failing test for schema creation**

Create `tests/conftest.py`:

```python
import sqlite3
import pytest


@pytest.fixture
def db_conn():
    """In-memory SQLite database for testing."""
    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row
    yield conn
    conn.close()
```

Create `tests/test_db.py`:

```python
from scripts.db import init_db, get_connection


def test_init_db_creates_all_tables(db_conn):
    init_db(db_conn)
    cursor = db_conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    )
    tables = [row["name"] for row in cursor.fetchall()]
    assert "drugs" in tables
    assert "drug_targets" in tables
    assert "keloid_targets" in tables


def test_init_db_keloid_targets_has_ensembl_id(db_conn):
    init_db(db_conn)
    cursor = db_conn.execute("PRAGMA table_info(keloid_targets)")
    columns = [row["name"] for row in cursor.fetchall()]
    assert "ensembl_id" in columns
    assert "gene_symbol" in columns


def test_init_db_gene_symbol_is_unique(db_conn):
    init_db(db_conn)
    db_conn.execute(
        "INSERT INTO keloid_targets (gene_symbol, target_name, pathway, evidence_type, evidence_strength, source)"
        " VALUES ('MTOR', 'mTOR', 'PI3K/AKT/mTOR', 'functional', 0.7, 'test')"
    )
    import sqlite3 as _sqlite3
    with pytest.raises(_sqlite3.IntegrityError):
        db_conn.execute(
            "INSERT INTO keloid_targets (gene_symbol, target_name, pathway, evidence_type, evidence_strength, source)"
            " VALUES ('MTOR', 'mTOR duplicate', 'PI3K/AKT/mTOR', 'functional', 0.5, 'test')"
        )


def test_get_connection_creates_file_db(tmp_path):
    db_path = tmp_path / "test.db"
    conn = get_connection(str(db_path))
    init_db(conn)
    cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
    assert len(cursor.fetchall()) == 3
    conn.close()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd ~/Coding\ Projects/Keloid
python3 -m pytest tests/test_db.py -v
```

Expected: `ModuleNotFoundError: No module named 'scripts'` or `ImportError`

- [ ] **Step 3: Implement db.py**

Create `scripts/__init__.py` (empty file) and `scripts/db.py`:

```python
import sqlite3
import os

DB_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "keloid.db")


def get_connection(db_path=None):
    """Return a SQLite connection with Row factory enabled."""
    path = db_path or DB_PATH
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def init_db(conn):
    """Create all tables if they don't exist."""
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS keloid_targets (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            target_name TEXT NOT NULL,
            gene_symbol TEXT NOT NULL UNIQUE,
            ensembl_id TEXT,
            pathway TEXT NOT NULL,
            evidence_type TEXT NOT NULL,
            evidence_strength REAL NOT NULL,
            source TEXT NOT NULL,
            created_at TEXT DEFAULT (datetime('now'))
        );

        CREATE TABLE IF NOT EXISTS drugs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            drug_name TEXT NOT NULL,
            generic_name TEXT NOT NULL UNIQUE,
            approval_status TEXT NOT NULL DEFAULT 'approved',
            original_indication TEXT,
            mechanism_of_action TEXT,
            created_at TEXT DEFAULT (datetime('now'))
        );

        CREATE TABLE IF NOT EXISTS drug_targets (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            drug_id INTEGER NOT NULL REFERENCES drugs(id),
            gene_symbol TEXT NOT NULL,
            action_type TEXT,
            source TEXT NOT NULL,
            created_at TEXT DEFAULT (datetime('now'))
        );
    """)
    conn.commit()
```

Also create empty `scripts/__init__.py` and `tests/__init__.py`:

```bash
touch scripts/__init__.py tests/__init__.py
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
python3 -m pytest tests/test_db.py -v
```

Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add scripts/db.py scripts/__init__.py tests/
git commit -m "feat: database layer with schema init and connection helper"
```

---

## Task 2: Keloid Target Seeding (`scripts/01_seed_keloid_targets.py`)

**Files:**
- Create: `scripts/01_seed_keloid_targets.py`
- Create: `tests/test_seed_targets.py`
- Create: `tests/fixtures/opentargets_keloid.json`

This script does three things:
1. Load manually curated targets from `data/seed_targets.csv`
2. Fetch keloid-associated targets from OpenTargets GraphQL API
3. Merge both sets into `keloid_targets` table (CSV seeds take priority on conflict)

### Data flow:
```
seed_targets.csv ──→ parse_csv() ──→ [{gene_symbol, target_name, ...}]
                                            │
OpenTargets API ──→ fetch_opentargets() ──→ [{gene_symbol, ensembl_id, ...}]
                                            │
                                      merge_targets()
                                            │
                                            ▼
                                   INSERT INTO keloid_targets
                                   (INSERT OR IGNORE — CSV loaded first,
                                    so CSV seeds win on gene_symbol conflict)
```

- [ ] **Step 1: Create the test fixture**

Create `tests/fixtures/opentargets_keloid.json`:

```json
{
  "data": {
    "disease": {
      "id": "EFO_0004212",
      "name": "Keloid",
      "associatedTargets": {
        "count": 5,
        "rows": [
          {
            "target": {
              "id": "ENSG00000174307",
              "approvedSymbol": "PHLDA3",
              "approvedName": "pleckstrin homology like domain family A member 3"
            },
            "score": 0.52,
            "datatypeScores": [
              {"id": "genetic_association", "score": 0.86}
            ]
          },
          {
            "target": {
              "id": "ENSG00000069869",
              "approvedSymbol": "NEDD4",
              "approvedName": "NEDD4 E3 ubiquitin protein ligase"
            },
            "score": 0.49,
            "datatypeScores": [
              {"id": "literature", "score": 0.54},
              {"id": "genetic_association", "score": 0.78}
            ]
          },
          {
            "target": {
              "id": "ENSG00000105329",
              "approvedSymbol": "TGFB1",
              "approvedName": "transforming growth factor beta 1"
            },
            "score": 0.35,
            "datatypeScores": [
              {"id": "literature", "score": 0.70}
            ]
          }
        ]
      }
    }
  }
}
```

Note: TGFB1 appears in BOTH the fixture AND the seed CSV. This tests the merge/dedup behavior.

- [ ] **Step 2: Write failing tests**

Create `tests/test_seed_targets.py`:

```python
import csv
import io
import json
import os
import pytest
from scripts.db import init_db
from scripts.seed_keloid_targets import (
    parse_seed_csv,
    parse_opentargets_response,
    insert_targets,
)


@pytest.fixture
def sample_csv():
    return """gene_symbol,target_name,pathway,evidence_type,evidence_strength,source
TGFB1,transforming growth factor beta 1,TGF-β/Smad,functional,0.7,manual_curation
MTOR,mechanistic target of rapamycin kinase,PI3K/AKT/mTOR,functional,0.7,manual_curation"""


@pytest.fixture
def opentargets_response():
    fixture_path = os.path.join(
        os.path.dirname(__file__), "fixtures", "opentargets_keloid.json"
    )
    with open(fixture_path) as f:
        return json.load(f)


class TestParseSeedCsv:
    def test_returns_list_of_dicts(self, sample_csv):
        targets = parse_seed_csv(io.StringIO(sample_csv))
        assert len(targets) == 2
        assert targets[0]["gene_symbol"] == "TGFB1"
        assert targets[1]["gene_symbol"] == "MTOR"

    def test_evidence_strength_is_float(self, sample_csv):
        targets = parse_seed_csv(io.StringIO(sample_csv))
        assert isinstance(targets[0]["evidence_strength"], float)
        assert targets[0]["evidence_strength"] == 0.7

    def test_empty_csv_returns_empty_list(self):
        empty = "gene_symbol,target_name,pathway,evidence_type,evidence_strength,source\n"
        targets = parse_seed_csv(io.StringIO(empty))
        assert targets == []


class TestParseOpentargetsResponse:
    def test_extracts_gene_symbols(self, opentargets_response):
        targets = parse_opentargets_response(opentargets_response)
        symbols = [t["gene_symbol"] for t in targets]
        assert "PHLDA3" in symbols
        assert "NEDD4" in symbols
        assert "TGFB1" in symbols

    def test_includes_ensembl_id(self, opentargets_response):
        targets = parse_opentargets_response(opentargets_response)
        phlda3 = [t for t in targets if t["gene_symbol"] == "PHLDA3"][0]
        assert phlda3["ensembl_id"] == "ENSG00000174307"

    def test_uses_opentargets_score(self, opentargets_response):
        targets = parse_opentargets_response(opentargets_response)
        phlda3 = [t for t in targets if t["gene_symbol"] == "PHLDA3"][0]
        assert phlda3["evidence_strength"] == 0.52

    def test_empty_response_returns_empty_list(self):
        empty = {"data": {"disease": None}}
        targets = parse_opentargets_response(empty)
        assert targets == []

    def test_zero_targets_returns_empty_list(self):
        zero = {"data": {"disease": {"associatedTargets": {"count": 0, "rows": []}}}}
        targets = parse_opentargets_response(zero)
        assert targets == []


class TestInsertTargets:
    def test_inserts_csv_targets(self, db_conn, sample_csv):
        init_db(db_conn)
        csv_targets = parse_seed_csv(io.StringIO(sample_csv))
        insert_targets(db_conn, csv_targets, [])
        cursor = db_conn.execute("SELECT COUNT(*) as cnt FROM keloid_targets")
        assert cursor.fetchone()["cnt"] == 2

    def test_csv_seeds_win_on_conflict(self, db_conn, sample_csv, opentargets_response):
        """CSV seeds are loaded first. When OpenTargets has the same gene_symbol
        (TGFB1), the CSV version should be kept (evidence_strength=0.7, not 0.35)."""
        init_db(db_conn)
        csv_targets = parse_seed_csv(io.StringIO(sample_csv))
        api_targets = parse_opentargets_response(opentargets_response)
        insert_targets(db_conn, csv_targets, api_targets)

        cursor = db_conn.execute(
            "SELECT evidence_strength, source FROM keloid_targets WHERE gene_symbol = 'TGFB1'"
        )
        row = cursor.fetchone()
        assert row["evidence_strength"] == 0.7
        assert row["source"] == "manual_curation"

    def test_api_targets_added_when_no_conflict(self, db_conn, sample_csv, opentargets_response):
        init_db(db_conn)
        csv_targets = parse_seed_csv(io.StringIO(sample_csv))
        api_targets = parse_opentargets_response(opentargets_response)
        insert_targets(db_conn, csv_targets, api_targets)

        cursor = db_conn.execute(
            "SELECT gene_symbol FROM keloid_targets ORDER BY gene_symbol"
        )
        symbols = [row["gene_symbol"] for row in cursor.fetchall()]
        # MTOR (csv only), NEDD4 (api only), PHLDA3 (api only), TGFB1 (both — csv wins)
        assert symbols == ["MTOR", "NEDD4", "PHLDA3", "TGFB1"]

    def test_total_count_with_overlap(self, db_conn, sample_csv, opentargets_response):
        init_db(db_conn)
        csv_targets = parse_seed_csv(io.StringIO(sample_csv))
        api_targets = parse_opentargets_response(opentargets_response)
        insert_targets(db_conn, csv_targets, api_targets)
        cursor = db_conn.execute("SELECT COUNT(*) as cnt FROM keloid_targets")
        # 2 from CSV + 2 new from API (TGFB1 is a duplicate) = 4
        assert cursor.fetchone()["cnt"] == 4
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
python3 -m pytest tests/test_seed_targets.py -v
```

Expected: `ImportError: cannot import name 'parse_seed_csv' from 'scripts.seed_keloid_targets'`

- [ ] **Step 4: Implement 01_seed_keloid_targets.py**

Create `scripts/01_seed_keloid_targets.py`:

```python
"""Seed keloid_targets table from curated CSV + OpenTargets API."""

import csv
import json
import logging
import os
import requests
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
KELOID_DISEASE_ID = "EFO_0004212"
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
SEED_CSV = os.path.join(DATA_DIR, "seed_targets.csv")

# Fetch up to 500 targets per page. OpenTargets has ~1600 for keloid
# but most have very low scores. We take the top 500.
PAGE_SIZE = 500

OPENTARGETS_QUERY = """
query keloidTargets($diseaseId: String!, $size: Int!) {
  disease(efoId: $diseaseId) {
    id
    name
    associatedTargets(page: {index: 0, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""


def parse_seed_csv(file_obj):
    """Parse the seed CSV into a list of target dicts."""
    reader = csv.DictReader(file_obj)
    targets = []
    for row in reader:
        targets.append({
            "gene_symbol": row["gene_symbol"],
            "target_name": row["target_name"],
            "ensembl_id": None,
            "pathway": row["pathway"],
            "evidence_type": row["evidence_type"],
            "evidence_strength": float(row["evidence_strength"]),
            "source": row["source"],
        })
    return targets


def parse_opentargets_response(response_json):
    """Parse OpenTargets GraphQL response into a list of target dicts."""
    disease = response_json.get("data", {}).get("disease")
    if disease is None:
        logger.warning("OpenTargets returned no disease data for keloid query")
        return []

    assoc = disease.get("associatedTargets", {})
    rows = assoc.get("rows", [])
    if not rows:
        logger.warning("OpenTargets returned 0 associated targets for keloid")
        return []

    targets = []
    for row in rows:
        target = row["target"]
        targets.append({
            "gene_symbol": target["approvedSymbol"],
            "target_name": target["approvedName"],
            "ensembl_id": target["id"],
            "pathway": "OpenTargets_association",
            "evidence_type": "computed",
            "evidence_strength": round(row["score"], 4),
            "source": "opentargets",
        })
    return targets


def fetch_opentargets(cache_path=None):
    """Fetch keloid targets from OpenTargets. Uses cache if available."""
    cache = cache_path or os.path.join(RAW_DIR, "opentargets_keloid.json")

    if os.path.exists(cache):
        logger.info(f"Loading cached OpenTargets response from {cache}")
        with open(cache) as f:
            return json.load(f)

    logger.info("Fetching keloid targets from OpenTargets API...")
    resp = requests.post(
        OPENTARGETS_URL,
        json={
            "query": OPENTARGETS_QUERY,
            "variables": {"diseaseId": KELOID_DISEASE_ID, "size": PAGE_SIZE},
        },
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()

    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Cached OpenTargets response to {cache}")

    return data


def insert_targets(conn, csv_targets, api_targets):
    """Insert targets into keloid_targets. CSV seeds are inserted first
    so they win on gene_symbol UNIQUE conflict (INSERT OR IGNORE)."""
    sql = """
        INSERT OR IGNORE INTO keloid_targets
        (gene_symbol, target_name, ensembl_id, pathway, evidence_type, evidence_strength, source)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    """
    for t in csv_targets:
        conn.execute(sql, (
            t["gene_symbol"], t["target_name"], t["ensembl_id"],
            t["pathway"], t["evidence_type"], t["evidence_strength"], t["source"],
        ))

    for t in api_targets:
        conn.execute(sql, (
            t["gene_symbol"], t["target_name"], t["ensembl_id"],
            t["pathway"], t["evidence_type"], t["evidence_strength"], t["source"],
        ))

    conn.commit()
    cursor = conn.execute("SELECT COUNT(*) as cnt FROM keloid_targets")
    count = cursor.fetchone()["cnt"]
    logger.info(f"Inserted {count} keloid targets")


def main():
    conn = get_connection()
    init_db(conn)

    # Clear existing targets for idempotent re-runs
    conn.execute("DELETE FROM keloid_targets")
    conn.commit()

    # Load CSV seeds
    with open(SEED_CSV) as f:
        csv_targets = parse_seed_csv(f)
    logger.info(f"Loaded {len(csv_targets)} targets from seed CSV")

    # Fetch from OpenTargets
    api_response = fetch_opentargets()
    api_targets = parse_opentargets_response(api_response)
    logger.info(f"Parsed {len(api_targets)} targets from OpenTargets")

    # Merge and insert
    insert_targets(conn, csv_targets, api_targets)
    conn.close()


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
python3 -m pytest tests/test_seed_targets.py -v
```

Expected: 10 passed

- [ ] **Step 6: Commit**

```bash
git add scripts/01_seed_keloid_targets.py tests/test_seed_targets.py tests/fixtures/opentargets_keloid.json
git commit -m "feat: keloid target seeding from CSV + OpenTargets API"
```

---

## Task 3: Drug Ingestion (`scripts/02_ingest_drugs.py`)

**Files:**
- Create: `scripts/02_ingest_drugs.py`
- Create: `tests/test_ingest_drugs.py`
- Create: `tests/fixtures/dgidb_interactions.json`

This script queries DGIdb for drug-gene interactions against every gene in `keloid_targets`, filters to FDA-approved drugs, and inserts into `drugs` + `drug_targets`.

### Data flow:
```
keloid_targets (SQLite) ──→ list of gene_symbols
                                    │
                                    ▼
DGIdb GraphQL API ──→ fetch per gene ──→ raw JSON cache
                                    │
                                    ▼
                           parse_dgidb_response()
                                    │
                           Filter: drug.approved == true
                                    │
                                    ▼
                           Dedup drugs by generic_name (lowercased)
                                    │
                              ┌─────┴─────┐
                              ▼           ▼
                     INSERT drugs    INSERT drug_targets
```

- [ ] **Step 1: Create the test fixture**

Create `tests/fixtures/dgidb_interactions.json`:

```json
{
  "data": {
    "genes": {
      "nodes": [
        {
          "name": "MTOR",
          "longName": "mechanistic target of rapamycin kinase",
          "interactions": [
            {
              "drug": {
                "name": "EVEROLIMUS",
                "approved": true,
                "conceptId": "rxcui:141704"
              },
              "interactionScore": 0.21,
              "interactionTypes": [
                {"type": "inhibitor", "directionality": "INHIBITORY"}
              ]
            },
            {
              "drug": {
                "name": "SIROLIMUS",
                "approved": true,
                "conceptId": "rxcui:35302"
              },
              "interactionScore": 0.09,
              "interactionTypes": [
                {"type": "inhibitor", "directionality": "INHIBITORY"}
              ]
            },
            {
              "drug": {
                "name": "VISTUSERTIB",
                "approved": false,
                "conceptId": "ncit:C88329"
              },
              "interactionScore": 0.09,
              "interactionTypes": [
                {"type": "inhibitor", "directionality": "INHIBITORY"}
              ]
            }
          ]
        },
        {
          "name": "EGFR",
          "longName": "epidermal growth factor receptor",
          "interactions": [
            {
              "drug": {
                "name": "ERLOTINIB",
                "approved": true,
                "conceptId": "rxcui:337535"
              },
              "interactionScore": 0.15,
              "interactionTypes": [
                {"type": "inhibitor", "directionality": "INHIBITORY"}
              ]
            },
            {
              "drug": {
                "name": "EVEROLIMUS",
                "approved": true,
                "conceptId": "rxcui:141704"
              },
              "interactionScore": 0.05,
              "interactionTypes": []
            }
          ]
        }
      ]
    }
  }
}
```

Note: EVEROLIMUS appears for both MTOR and EGFR. This tests drug deduplication and multi-target tracking.

- [ ] **Step 2: Write failing tests**

Create `tests/test_ingest_drugs.py`:

```python
import json
import os
import pytest
from scripts.db import init_db
from scripts.ingest_drugs import (
    parse_dgidb_response,
    insert_drugs_and_targets,
)


@pytest.fixture
def dgidb_response():
    fixture_path = os.path.join(
        os.path.dirname(__file__), "fixtures", "dgidb_interactions.json"
    )
    with open(fixture_path) as f:
        return json.load(f)


class TestParseDgidbResponse:
    def test_filters_to_approved_only(self, dgidb_response):
        drugs, interactions = parse_dgidb_response(dgidb_response)
        drug_names = [d["drug_name"] for d in drugs]
        assert "VISTUSERTIB" not in drug_names
        assert "EVEROLIMUS" in drug_names
        assert "SIROLIMUS" in drug_names

    def test_deduplicates_drugs(self, dgidb_response):
        """EVEROLIMUS appears for both MTOR and EGFR. Should appear once in drugs."""
        drugs, interactions = parse_dgidb_response(dgidb_response)
        everolimus_count = sum(1 for d in drugs if d["drug_name"] == "EVEROLIMUS")
        assert everolimus_count == 1

    def test_preserves_all_interactions(self, dgidb_response):
        """EVEROLIMUS should have interactions with both MTOR and EGFR."""
        drugs, interactions = parse_dgidb_response(dgidb_response)
        everolimus_genes = [
            i["gene_symbol"] for i in interactions if i["drug_name"] == "EVEROLIMUS"
        ]
        assert "MTOR" in everolimus_genes
        assert "EGFR" in everolimus_genes

    def test_extracts_action_type(self, dgidb_response):
        drugs, interactions = parse_dgidb_response(dgidb_response)
        sirolimus_int = [i for i in interactions if i["drug_name"] == "SIROLIMUS"][0]
        assert sirolimus_int["action_type"] == "inhibitor"

    def test_empty_action_type_defaults_to_unknown(self, dgidb_response):
        """EVEROLIMUS + EGFR has no interactionTypes. Should default to 'unknown'."""
        drugs, interactions = parse_dgidb_response(dgidb_response)
        ev_egfr = [
            i for i in interactions
            if i["drug_name"] == "EVEROLIMUS" and i["gene_symbol"] == "EGFR"
        ][0]
        assert ev_egfr["action_type"] == "unknown"

    def test_empty_response_returns_empty(self):
        empty = {"data": {"genes": {"nodes": []}}}
        drugs, interactions = parse_dgidb_response(empty)
        assert drugs == []
        assert interactions == []


class TestInsertDrugsAndTargets:
    def test_inserts_drugs(self, db_conn, dgidb_response):
        init_db(db_conn)
        drugs, interactions = parse_dgidb_response(dgidb_response)
        insert_drugs_and_targets(db_conn, drugs, interactions)
        cursor = db_conn.execute("SELECT COUNT(*) as cnt FROM drugs")
        # EVEROLIMUS, SIROLIMUS, ERLOTINIB = 3
        assert cursor.fetchone()["cnt"] == 3

    def test_inserts_drug_targets(self, db_conn, dgidb_response):
        init_db(db_conn)
        drugs, interactions = parse_dgidb_response(dgidb_response)
        insert_drugs_and_targets(db_conn, drugs, interactions)
        cursor = db_conn.execute("SELECT COUNT(*) as cnt FROM drug_targets")
        # EVEROLIMUS→MTOR, EVEROLIMUS→EGFR, SIROLIMUS→MTOR, ERLOTINIB→EGFR = 4
        assert cursor.fetchone()["cnt"] == 4

    def test_drug_target_references_correct_drug_id(self, db_conn, dgidb_response):
        init_db(db_conn)
        drugs, interactions = parse_dgidb_response(dgidb_response)
        insert_drugs_and_targets(db_conn, drugs, interactions)
        cursor = db_conn.execute("""
            SELECT d.drug_name, dt.gene_symbol
            FROM drug_targets dt JOIN drugs d ON dt.drug_id = d.id
            WHERE d.drug_name = 'EVEROLIMUS'
            ORDER BY dt.gene_symbol
        """)
        rows = cursor.fetchall()
        genes = [r["gene_symbol"] for r in rows]
        assert genes == ["EGFR", "MTOR"]
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
python3 -m pytest tests/test_ingest_drugs.py -v
```

Expected: `ImportError`

- [ ] **Step 4: Implement 02_ingest_drugs.py**

Create `scripts/02_ingest_drugs.py`:

```python
"""Ingest drug-gene interactions from DGIdb for all keloid target genes."""

import json
import logging
import os
import requests
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DGIDB_URL = "https://dgidb.org/api/graphql"
RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw")

DGIDB_QUERY = """
query drugInteractions($genes: [String!]!) {
  genes(names: $genes) {
    nodes {
      name
      longName
      interactions {
        drug {
          name
          approved
          conceptId
        }
        interactionScore
        interactionTypes {
          type
          directionality
        }
      }
    }
  }
}
"""


def fetch_dgidb(gene_symbols, cache_path=None):
    """Fetch drug-gene interactions from DGIdb. Uses cache if available."""
    cache = cache_path or os.path.join(RAW_DIR, "dgidb_interactions.json")

    if os.path.exists(cache):
        logger.info(f"Loading cached DGIdb response from {cache}")
        with open(cache) as f:
            return json.load(f)

    logger.info(f"Fetching interactions for {len(gene_symbols)} genes from DGIdb...")
    resp = requests.post(
        DGIDB_URL,
        json={"query": DGIDB_QUERY, "variables": {"genes": gene_symbols}},
        timeout=60,
    )
    resp.raise_for_status()
    data = resp.json()

    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Cached DGIdb response to {cache}")

    return data


def parse_dgidb_response(response_json):
    """Parse DGIdb response into deduplicated drugs list and interactions list.
    Filters to FDA-approved drugs only."""
    nodes = response_json.get("data", {}).get("genes", {}).get("nodes", [])

    seen_drugs = {}  # lowercase drug_name → drug dict
    interactions = []

    for gene_node in nodes:
        gene_symbol = gene_node["name"]
        for interaction in gene_node.get("interactions", []):
            drug = interaction["drug"]
            if not drug.get("approved", False):
                continue

            drug_name = drug["name"]
            drug_key = drug_name.lower()

            if drug_key not in seen_drugs:
                seen_drugs[drug_key] = {
                    "drug_name": drug_name,
                    "generic_name": drug_name.lower(),
                    "approval_status": "approved",
                    "original_indication": None,
                    "mechanism_of_action": None,
                }

            action_types = interaction.get("interactionTypes", [])
            action_type = action_types[0]["type"] if action_types else "unknown"

            interactions.append({
                "drug_name": drug_name,
                "gene_symbol": gene_symbol,
                "action_type": action_type,
                "source": "DGIdb",
            })

    drugs = list(seen_drugs.values())
    return drugs, interactions


def insert_drugs_and_targets(conn, drugs, interactions):
    """Insert drugs and drug_targets into SQLite.
    Drugs are inserted first, then interactions reference drug IDs."""
    drug_id_map = {}  # drug_name → id
    for d in drugs:
        cursor = conn.execute(
            """INSERT OR IGNORE INTO drugs (drug_name, generic_name, approval_status,
               original_indication, mechanism_of_action)
               VALUES (?, ?, ?, ?, ?)""",
            (d["drug_name"], d["generic_name"], d["approval_status"],
             d["original_indication"], d["mechanism_of_action"]),
        )
        if cursor.lastrowid:
            drug_id_map[d["drug_name"]] = cursor.lastrowid
        else:
            # Already existed — fetch its ID
            row = conn.execute(
                "SELECT id FROM drugs WHERE generic_name = ?", (d["generic_name"],)
            ).fetchone()
            drug_id_map[d["drug_name"]] = row["id"]

    for i in interactions:
        drug_id = drug_id_map[i["drug_name"]]
        conn.execute(
            """INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source)
               VALUES (?, ?, ?, ?)""",
            (drug_id, i["gene_symbol"], i["action_type"], i["source"]),
        )

    conn.commit()
    drug_count = conn.execute("SELECT COUNT(*) as cnt FROM drugs").fetchone()["cnt"]
    target_count = conn.execute("SELECT COUNT(*) as cnt FROM drug_targets").fetchone()["cnt"]
    logger.info(f"Inserted {drug_count} drugs with {target_count} drug-target interactions")


def main():
    conn = get_connection()
    init_db(conn)

    # Clear existing drug data for idempotent re-runs
    conn.execute("DELETE FROM drug_targets")
    conn.execute("DELETE FROM drugs")
    conn.commit()

    # Get all gene symbols from keloid_targets
    cursor = conn.execute("SELECT gene_symbol FROM keloid_targets")
    gene_symbols = [row["gene_symbol"] for row in cursor.fetchall()]

    if not gene_symbols:
        logger.warning("No keloid targets found. Run 01_seed_keloid_targets.py first.")
        conn.close()
        return

    logger.info(f"Querying DGIdb for {len(gene_symbols)} keloid target genes")

    # Fetch from DGIdb
    response = fetch_dgidb(gene_symbols)
    drugs, interactions = parse_dgidb_response(response)
    logger.info(f"Found {len(drugs)} approved drugs with {len(interactions)} interactions")

    # Insert
    insert_drugs_and_targets(conn, drugs, interactions)
    conn.close()


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
python3 -m pytest tests/test_ingest_drugs.py -v
```

Expected: 10 passed

- [ ] **Step 6: Commit**

```bash
git add scripts/02_ingest_drugs.py tests/test_ingest_drugs.py tests/fixtures/dgidb_interactions.json
git commit -m "feat: drug ingestion from DGIdb with approval filtering and dedup"
```

---

## Task 4: Scoring Engine (`scripts/03_score_candidates.py`)

**Files:**
- Create: `scripts/03_score_candidates.py`
- Create: `tests/test_scoring.py`

This is the core of the pipeline. It joins `drug_targets` against `keloid_targets`, counts pathway overlaps per drug, and computes a composite score.

### Scoring formula:
```
v1_score = (
    normalized_overlap   × 0.40 +     # pathway_overlap_count / max_overlap_count
    avg_evidence_strength × 0.40 +     # average across matched targets
    multi_target_bonus    × 0.20       # 1.0 if 3+ distinct pathways, else 0.0
)
```

Note: `pathway_overlap_count` is normalized to 0.0–1.0 by dividing by the maximum overlap count across all drugs. This keeps the score components on the same scale.

- [ ] **Step 1: Write failing tests**

Create `tests/test_scoring.py`:

```python
import pytest
from scripts.db import init_db
from scripts.score_candidates import (
    query_overlaps,
    compute_scores,
    normalize_overlap_counts,
)


@pytest.fixture
def populated_db(db_conn):
    """DB with keloid targets and drugs that have known overlaps."""
    init_db(db_conn)

    # Insert keloid targets across 3 pathways
    db_conn.executemany(
        """INSERT INTO keloid_targets (gene_symbol, target_name, pathway, evidence_type, evidence_strength, source)
           VALUES (?, ?, ?, ?, ?, ?)""",
        [
            ("MTOR", "mTOR kinase", "PI3K/AKT/mTOR", "functional", 0.7, "manual"),
            ("PIK3CA", "PI3K alpha", "PI3K/AKT/mTOR", "functional", 0.6, "manual"),
            ("EGFR", "EGF receptor", "EGFR", "expression", 0.8, "manual"),
            ("TGFB1", "TGF-beta 1", "TGF-β/Smad", "functional", 0.9, "manual"),
            ("JAK2", "Janus kinase 2", "JAK/STAT", "functional", 0.5, "opentargets"),
        ],
    )

    # Drug A: hits 3 pathways (MTOR, EGFR, TGFB1) — should get multi-target bonus
    db_conn.execute(
        "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (1, 'DRUG_A', 'drug_a', 'approved')"
    )
    db_conn.executemany(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (?, ?, ?, ?)",
        [
            (1, "MTOR", "inhibitor", "DGIdb"),
            (1, "EGFR", "inhibitor", "DGIdb"),
            (1, "TGFB1", "inhibitor", "DGIdb"),
        ],
    )

    # Drug B: hits 1 pathway (MTOR only) — no multi-target bonus
    db_conn.execute(
        "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (2, 'DRUG_B', 'drug_b', 'approved')"
    )
    db_conn.execute(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (2, 'MTOR', 'inhibitor', 'DGIdb')"
    )

    # Drug C: hits 2 genes in SAME pathway (MTOR, PIK3CA both in PI3K/AKT/mTOR) — 1 distinct pathway
    db_conn.execute(
        "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (3, 'DRUG_C', 'drug_c', 'approved')"
    )
    db_conn.executemany(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (?, ?, ?, ?)",
        [
            (3, "MTOR", "inhibitor", "DGIdb"),
            (3, "PIK3CA", "inhibitor", "DGIdb"),
        ],
    )

    # Drug D: targets JAK2 only (low evidence pathway)
    db_conn.execute(
        "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (4, 'DRUG_D', 'drug_d', 'approved')"
    )
    db_conn.execute(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (4, 'JAK2', 'inhibitor', 'DGIdb')"
    )

    db_conn.commit()
    return db_conn


class TestQueryOverlaps:
    def test_returns_all_matching_drugs(self, populated_db):
        overlaps = query_overlaps(populated_db)
        drug_names = set(o["drug_name"] for o in overlaps)
        assert drug_names == {"DRUG_A", "DRUG_B", "DRUG_C", "DRUG_D"}

    def test_drug_a_has_three_overlaps(self, populated_db):
        overlaps = query_overlaps(populated_db)
        drug_a = [o for o in overlaps if o["drug_name"] == "DRUG_A"]
        assert len(drug_a) == 3

    def test_no_overlaps_when_gene_not_in_keloid_targets(self, populated_db):
        """Drug targeting a gene NOT in keloid_targets should not appear."""
        populated_db.execute(
            "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (99, 'ORPHAN', 'orphan', 'approved')"
        )
        populated_db.execute(
            "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (99, 'NONEXISTENT', 'inhibitor', 'DGIdb')"
        )
        populated_db.commit()
        overlaps = query_overlaps(populated_db)
        orphan = [o for o in overlaps if o["drug_name"] == "ORPHAN"]
        assert orphan == []


class TestComputeScores:
    def test_drug_a_gets_multi_target_bonus(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_a = [s for s in scores if s["drug_name"] == "DRUG_A"][0]
        assert drug_a["multi_target_bonus"] == 1.0

    def test_drug_b_no_multi_target_bonus(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_b = [s for s in scores if s["drug_name"] == "DRUG_B"][0]
        assert drug_b["multi_target_bonus"] == 0.0

    def test_drug_c_counts_one_distinct_pathway(self, populated_db):
        """MTOR and PIK3CA are both PI3K/AKT/mTOR — should count as 1 pathway."""
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_c = [s for s in scores if s["drug_name"] == "DRUG_C"][0]
        assert drug_c["pathway_overlap_count"] == 1

    def test_drug_a_ranked_above_drug_b(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        scores.sort(key=lambda s: s["composite_score"], reverse=True)
        names = [s["drug_name"] for s in scores]
        assert names.index("DRUG_A") < names.index("DRUG_B")

    def test_avg_evidence_strength(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_b = [s for s in scores if s["drug_name"] == "DRUG_B"][0]
        # DRUG_B only hits MTOR (evidence_strength=0.7)
        assert drug_b["avg_evidence_strength"] == 0.7

    def test_drug_d_low_evidence(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_d = [s for s in scores if s["drug_name"] == "DRUG_D"][0]
        assert drug_d["avg_evidence_strength"] == 0.5


class TestNormalizeOverlapCounts:
    def test_max_overlap_normalized_to_one(self):
        scores = [
            {"drug_name": "A", "pathway_overlap_count": 4},
            {"drug_name": "B", "pathway_overlap_count": 2},
            {"drug_name": "C", "pathway_overlap_count": 1},
        ]
        normalized = normalize_overlap_counts(scores)
        assert normalized[0]["normalized_overlap"] == 1.0
        assert normalized[1]["normalized_overlap"] == 0.5
        assert normalized[2]["normalized_overlap"] == 0.25

    def test_single_drug_normalized_to_one(self):
        scores = [{"drug_name": "A", "pathway_overlap_count": 1}]
        normalized = normalize_overlap_counts(scores)
        assert normalized[0]["normalized_overlap"] == 1.0

    def test_empty_list(self):
        assert normalize_overlap_counts([]) == []
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
python3 -m pytest tests/test_scoring.py -v
```

Expected: `ImportError`

- [ ] **Step 3: Implement 03_score_candidates.py**

Create `scripts/03_score_candidates.py`:

```python
"""Score drug candidates by keloid pathway overlap."""

import logging
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Scoring weights (must sum to 1.0)
W_OVERLAP = 0.40
W_EVIDENCE = 0.40
W_MULTITARGET = 0.20

# Minimum distinct pathways to get multi-target bonus
MULTI_TARGET_THRESHOLD = 3


def query_overlaps(conn):
    """Find all drugs that target keloid-associated genes.
    Returns one row per drug-gene match with pathway and evidence info."""
    cursor = conn.execute("""
        SELECT d.id as drug_id, d.drug_name, d.generic_name,
               d.original_indication, d.mechanism_of_action,
               dt.gene_symbol, dt.action_type,
               kt.pathway, kt.evidence_strength
        FROM drugs d
        JOIN drug_targets dt ON d.id = dt.drug_id
        JOIN keloid_targets kt ON dt.gene_symbol = kt.gene_symbol
        WHERE d.approval_status = 'approved'
        ORDER BY d.drug_name, kt.pathway
    """)
    return [dict(row) for row in cursor.fetchall()]


def normalize_overlap_counts(scores):
    """Normalize pathway_overlap_count to 0.0–1.0 range."""
    if not scores:
        return []
    max_count = max(s["pathway_overlap_count"] for s in scores)
    for s in scores:
        s["normalized_overlap"] = s["pathway_overlap_count"] / max_count
    return scores


def compute_scores(overlaps):
    """Aggregate overlaps into per-drug scores.

    For each drug:
    - pathway_overlap_count: number of DISTINCT keloid pathways hit
    - avg_evidence_strength: mean evidence_strength across matched targets
    - multi_target_bonus: 1.0 if 3+ distinct pathways, else 0.0
    """
    # Group by drug
    drug_data = {}
    for row in overlaps:
        name = row["drug_name"]
        if name not in drug_data:
            drug_data[name] = {
                "drug_name": name,
                "generic_name": row["generic_name"],
                "original_indication": row["original_indication"],
                "mechanism_of_action": row["mechanism_of_action"],
                "pathways": set(),
                "genes": [],
                "evidence_values": [],
                "gene_details": [],
            }
        drug_data[name]["pathways"].add(row["pathway"])
        drug_data[name]["genes"].append(row["gene_symbol"])
        drug_data[name]["evidence_values"].append(row["evidence_strength"])
        drug_data[name]["gene_details"].append({
            "gene_symbol": row["gene_symbol"],
            "pathway": row["pathway"],
            "action_type": row["action_type"],
            "evidence_strength": row["evidence_strength"],
        })

    scores = []
    for name, data in drug_data.items():
        pathway_count = len(data["pathways"])
        avg_evidence = sum(data["evidence_values"]) / len(data["evidence_values"])
        multi_bonus = 1.0 if pathway_count >= MULTI_TARGET_THRESHOLD else 0.0

        scores.append({
            "drug_name": name,
            "generic_name": data["generic_name"],
            "original_indication": data["original_indication"],
            "mechanism_of_action": data["mechanism_of_action"],
            "pathway_overlap_count": pathway_count,
            "pathways": sorted(data["pathways"]),
            "avg_evidence_strength": round(avg_evidence, 4),
            "multi_target_bonus": multi_bonus,
            "gene_details": data["gene_details"],
            "composite_score": 0.0,  # computed after normalization
        })

    # Normalize and compute final score
    scores = normalize_overlap_counts(scores)
    for s in scores:
        s["composite_score"] = round(
            s["normalized_overlap"] * W_OVERLAP
            + s["avg_evidence_strength"] * W_EVIDENCE
            + s["multi_target_bonus"] * W_MULTITARGET,
            4,
        )

    scores.sort(key=lambda s: s["composite_score"], reverse=True)
    return scores


def main():
    conn = get_connection()

    overlaps = query_overlaps(conn)
    if not overlaps:
        logger.warning("No drug-keloid overlaps found. Check that steps 01 and 02 ran successfully.")
        conn.close()
        return

    scores = compute_scores(overlaps)
    logger.info(f"Scored {len(scores)} drug candidates")

    # Print top 10 summary to stdout
    print(f"\n{'Rank':<6}{'Drug':<25}{'Score':<8}{'Pathways':<10}{'Avg Evidence':<14}{'Multi-target'}")
    print("-" * 75)
    for i, s in enumerate(scores[:10], 1):
        print(
            f"{i:<6}{s['drug_name']:<25}{s['composite_score']:<8.4f}"
            f"{s['pathway_overlap_count']:<10}{s['avg_evidence_strength']:<14.4f}"
            f"{'YES' if s['multi_target_bonus'] else 'no'}"
        )

    # Save full scores as JSON for the report generator
    import json
    import os

    output_dir = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
    os.makedirs(output_dir, exist_ok=True)
    scores_path = os.path.join(output_dir, "scores.json")
    with open(scores_path, "w") as f:
        json.dump(scores, f, indent=2, default=list)
    logger.info(f"Full scores saved to {scores_path}")

    conn.close()


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
python3 -m pytest tests/test_scoring.py -v
```

Expected: 12 passed

- [ ] **Step 5: Commit**

```bash
git add scripts/03_score_candidates.py tests/test_scoring.py
git commit -m "feat: scoring engine with pathway overlap, evidence strength, and multi-target bonus"
```

---

## Task 5: Report Generator (`scripts/04_generate_report.py`)

**Files:**
- Create: `scripts/04_generate_report.py`

No tests for this task — it's pure Markdown templating with no branching logic worth testing (per plan decision #6).

- [ ] **Step 1: Implement 04_generate_report.py**

Create `scripts/04_generate_report.py`:

```python
"""Generate a Markdown report from scored drug candidates."""

import json
import logging
import os
from datetime import date

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

SCORES_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "scores.json")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output")


def generate_report(scores, top_n=20):
    """Generate Markdown report string from scored candidates."""
    today = date.today().isoformat()
    lines = [
        f"# Keloid Drug Repurposing Candidates",
        f"",
        f"*Generated {today} — {len(scores)} candidates scored*",
        f"",
        f"**Scoring formula:** `v1_score = overlap×0.40 + evidence×0.40 + multi_target×0.20`",
        f"",
        f"> This is a hypothesis generator, not medical advice. Discuss any candidate with a dermatologist before use.",
        f"",
        f"---",
        f"",
        f"## Top {min(top_n, len(scores))} Candidates",
        f"",
        f"| Rank | Drug | Score | Pathways Hit | Avg Evidence | Multi-target |",
        f"|------|------|-------|-------------|-------------|-------------|",
    ]

    for i, s in enumerate(scores[:top_n], 1):
        mt = "Yes" if s["multi_target_bonus"] else "-"
        lines.append(
            f"| {i} | {s['drug_name']} | {s['composite_score']:.4f} | "
            f"{s['pathway_overlap_count']} | {s['avg_evidence_strength']:.4f} | {mt} |"
        )

    lines.extend(["", "---", ""])

    # Detailed briefs for top 10
    lines.extend([f"## Detailed Briefs (Top 10)", ""])
    for i, s in enumerate(scores[:10], 1):
        lines.append(f"### {i}. {s['drug_name']}")
        lines.append("")
        if s.get("original_indication"):
            lines.append(f"**Approved for:** {s['original_indication']}")
        if s.get("mechanism_of_action"):
            lines.append(f"**Mechanism:** {s['mechanism_of_action']}")
        lines.append(f"**Composite score:** {s['composite_score']:.4f}")
        lines.append(f"**Distinct keloid pathways:** {s['pathway_overlap_count']} ({', '.join(s['pathways'])})")
        lines.append("")

        lines.append("| Gene | Pathway | Action | Evidence |")
        lines.append("|------|---------|--------|----------|")
        for g in s.get("gene_details", []):
            lines.append(
                f"| {g['gene_symbol']} | {g['pathway']} | "
                f"{g['action_type']} | {g['evidence_strength']:.2f} |"
            )
        lines.extend([
            "",
            "**Next steps:** Discuss with dermatologist. Search PubMed for "
            f"`\"{s['drug_name'].lower()}\" AND (keloid OR fibrosis)`. "
            "Check ClinicalTrials.gov for existing trials.",
            "",
            "---",
            "",
        ])

    return "\n".join(lines)


def main():
    if not os.path.exists(SCORES_PATH):
        logger.error(f"Scores file not found at {SCORES_PATH}. Run 03_score_candidates.py first.")
        return

    with open(SCORES_PATH) as f:
        scores = json.load(f)

    if not scores:
        logger.warning("No scored candidates found.")
        return

    report = generate_report(scores)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_path = os.path.join(OUTPUT_DIR, "report.md")
    with open(output_path, "w") as f:
        f.write(report)

    logger.info(f"Report written to {output_path} ({len(scores)} candidates)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add scripts/04_generate_report.py
git commit -m "feat: Markdown report generator with summary table and detailed briefs"
```

---

## Task 6: End-to-End Smoke Test

**Files:** None new — this is a manual verification step.

- [ ] **Step 1: Run all tests**

```bash
python3 -m pytest tests/ -v
```

Expected: All tests pass (approximately 32 tests across 4 test files).

- [ ] **Step 2: Run the full pipeline**

```bash
cd ~/Coding\ Projects/Keloid
./run_pipeline.sh
```

Expected output:
- Step 1 prints count of keloid targets loaded
- Step 2 prints count of drugs and interactions
- Step 3 prints top 10 ranked candidates table
- Step 4 confirms report written to `output/report.md`

- [ ] **Step 3: Inspect the report**

```bash
head -60 output/report.md
```

Verify: the report has a summary table and at least a few detailed briefs with drug names, pathways, and evidence.

- [ ] **Step 4: Commit any fixes, then final commit**

```bash
git add -A
git commit -m "chore: pipeline v1 complete — end-to-end verified"
```

---

## Self-Review Checklist

**Spec coverage:**
- Keloid pathway map from CSV + OpenTargets: Task 0 (CSV) + Task 2
- Drug-target ingestion from DGIdb: Task 3
- Pathway overlap scoring with multi-target bonus: Task 4
- Ranked Markdown report: Task 5
- Raw API caching: Built into Tasks 2 and 3 (fetch functions)
- Gene symbol normalization (Ensembl → HGNC): OpenTargets already returns `approvedSymbol` which IS the HGNC symbol. No separate normalization step needed — the API does it for us.
- Critical gap handling (empty API, unmapped IDs): Logging/warnings in Tasks 2 and 3

**Placeholder scan:** No TBD/TODO/placeholders found. All code blocks are complete.

**Type consistency:** `parse_seed_csv`, `parse_opentargets_response`, `parse_dgidb_response`, `query_overlaps`, `compute_scores`, `normalize_overlap_counts`, `generate_report` — all referenced consistently across tests and implementations.
