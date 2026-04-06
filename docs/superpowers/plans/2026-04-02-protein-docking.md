# Protein-Substrate Molecular Docking Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add structural validation to the keloid drug repurposing pipeline by docking top drug candidates into their target proteins using AutoDock Vina, filtering out false positives where a drug can't physically bind.

**Architecture:** New pipeline stage 06 fetches 3D structures for proteins (PDB/AlphaFold) and drugs (PubChem), prepares them for docking (PDBQT format via Open Babel), runs AutoDock Vina for each drug-target pair, and stores binding energy scores in SQLite. The existing composite score is updated to incorporate structural binding evidence.

**Tech Stack:** Python 3, AutoDock Vina (CLI binary), Open Babel (`brew install open-babel`), requests, SQLite. No conda or GPU required.

---

## Scope

**Proteins:** The 22 curated keloid target genes from `seed_targets.csv`. Of these, ~12 are kinases/enzymes with well-characterized binding pockets suitable for docking (EGFR, MTOR, PIK3CA, AKT1, JAK2, MAPK1, MAPK3, KDR, CYP24A1, STAT3, SMAD3, TNF). Structural proteins (COL1A1/A2, COL3A1) and secreted ligands (TGFB1/B2, WNT3A, VEGFA, IL6) are skipped — they lack drugable binding pockets.

**Drugs:** The 25 green/yellow severity tier drugs. Each drug is only docked against the curated target proteins it actually hits (from the `drug_targets` table).

**Estimated pairs:** ~60-80 drug-protein docking runs. At 2-10 min each on CPU, that's 2-13 hours total for the full batch.

---

## File Structure

```
scripts/
├── protein_structures.py    # Core: fetch PDB/AlphaFold structures, store in DB
├── drug_structures.py       # Core: fetch PubChem 3D SDF, convert to PDBQT
├── docking.py               # Core: run Vina docking, parse results
├── 06_structural_docking.py # Entry point: orchestrates fetch + dock + integrate
data/
├── structures/
│   ├── proteins/            # Downloaded .pdb files
│   ├── ligands/             # Downloaded .sdf and converted .pdbqt files
│   └── docking_results/     # Vina output poses
tests/
├── test_protein_structures.py
├── test_drug_structures.py
├── test_docking.py
├── fixtures/
│   ├── sample_pdb_search.json    # Cached PDB API response
│   ├── sample_alphafold.json     # Cached AlphaFold API response
│   ├── sample_pubchem_3d.json    # Cached PubChem conformer response
│   └── sample_vina_output.txt    # Sample Vina stdout
```

---

### Task 1: Install Dependencies and Verify

**Files:**
- Modify: `requirements.txt`

- [ ] **Step 1: Install Open Babel via Homebrew**

```bash
brew install open-babel
```

Verify: `obabel -V` should print version info.

- [ ] **Step 2: Download AutoDock Vina binary**

```bash
brew install autodock-vina
```

If not in Homebrew, download from https://github.com/ccsb-scripps/AutoDock-Vina/releases — get the `vina_1.2.5_mac_arm64` binary, place in `/usr/local/bin/vina` or add to PATH.

Verify: `vina --version` should print version.

- [ ] **Step 3: Install Python packages**

```bash
pip3 install meeko
```

Verify:

```bash
python3 -c "from meeko import MoleculePreparation; print('meeko OK')"
```

- [ ] **Step 4: Update requirements.txt**

Add `meeko>=0.5.0` to `requirements.txt`.

- [ ] **Step 5: Commit**

```bash
git add requirements.txt
git commit -m "chore: add meeko dependency for ligand preparation"
```

---

### Task 2: Extend Database Schema

**Files:**
- Modify: `scripts/db.py`
- Test: `tests/test_db.py`

- [ ] **Step 1: Write failing test for new tables**

Add to `tests/test_db.py`:

```python
from scripts.db import init_db


def test_protein_structures_table_created(db_conn):
    init_db(db_conn)
    cursor = db_conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='protein_structures'"
    )
    assert cursor.fetchone() is not None


def test_drug_structures_table_created(db_conn):
    init_db(db_conn)
    cursor = db_conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='drug_structures'"
    )
    assert cursor.fetchone() is not None


def test_docking_results_table_created(db_conn):
    init_db(db_conn)
    cursor = db_conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='docking_results'"
    )
    assert cursor.fetchone() is not None


def test_insert_protein_structure(db_conn):
    init_db(db_conn)
    db_conn.execute(
        """INSERT INTO protein_structures
           (gene_symbol, pdb_id, structure_source, resolution_angstroms, file_path)
           VALUES ('EGFR', '1M17', 'pdb_experimental', 2.6, 'data/structures/proteins/EGFR_1M17.pdb')"""
    )
    db_conn.commit()
    row = db_conn.execute("SELECT * FROM protein_structures WHERE gene_symbol='EGFR'").fetchone()
    assert row["pdb_id"] == "1M17"
    assert row["structure_source"] == "pdb_experimental"


def test_insert_docking_result(db_conn):
    init_db(db_conn)
    # Insert required foreign key rows
    db_conn.execute(
        """INSERT INTO protein_structures (id, gene_symbol, pdb_id, structure_source, file_path)
           VALUES (1, 'EGFR', '1M17', 'pdb_experimental', 'proteins/EGFR.pdb')"""
    )
    db_conn.execute(
        """INSERT INTO drug_structures (id, drug_name, pubchem_cid, smiles, sdf_path, pdbqt_path)
           VALUES (1, 'IMATINIB', '5291', 'CC1=CC=CC=C1', 'ligands/imatinib.sdf', 'ligands/imatinib.pdbqt')"""
    )
    db_conn.execute(
        """INSERT INTO docking_results
           (drug_structure_id, protein_structure_id, method, binding_energy_kcal)
           VALUES (1, 1, 'vina', -8.3)"""
    )
    db_conn.commit()
    row = db_conn.execute("SELECT * FROM docking_results").fetchone()
    assert row["binding_energy_kcal"] == -8.3
    assert row["method"] == "vina"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest tests/test_db.py -v`
Expected: FAIL — new tables don't exist yet.

- [ ] **Step 3: Add new tables to init_db**

In `scripts/db.py`, add these tables to the `init_db` function's `executescript` call, after the existing `drug_targets` table:

```python
        CREATE TABLE IF NOT EXISTS protein_structures (
            id INTEGER PRIMARY KEY,
            gene_symbol TEXT NOT NULL,
            pdb_id TEXT,
            alphafold_id TEXT,
            structure_source TEXT NOT NULL,
            resolution_angstroms REAL,
            binding_site_center_x REAL,
            binding_site_center_y REAL,
            binding_site_center_z REAL,
            binding_site_size REAL DEFAULT 22.0,
            file_path TEXT NOT NULL,
            created_at TEXT DEFAULT (datetime('now'))
        );

        CREATE TABLE IF NOT EXISTS drug_structures (
            id INTEGER PRIMARY KEY,
            drug_name TEXT NOT NULL,
            pubchem_cid TEXT,
            smiles TEXT,
            sdf_path TEXT,
            pdbqt_path TEXT,
            created_at TEXT DEFAULT (datetime('now'))
        );

        CREATE TABLE IF NOT EXISTS docking_results (
            id INTEGER PRIMARY KEY,
            drug_structure_id INTEGER NOT NULL REFERENCES drug_structures(id),
            protein_structure_id INTEGER NOT NULL REFERENCES protein_structures(id),
            method TEXT NOT NULL DEFAULT 'vina',
            binding_energy_kcal REAL,
            num_poses INTEGER,
            best_pose_path TEXT,
            created_at TEXT DEFAULT (datetime('now'))
        );
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m pytest tests/test_db.py -v`
Expected: All tests PASS.

- [ ] **Step 5: Run full test suite to check no regressions**

Run: `python3 -m pytest tests/ -v`
Expected: All 57 existing tests + new tests PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/db.py tests/test_db.py
git commit -m "feat: add protein_structures, drug_structures, docking_results tables"
```

---

### Task 3: Protein Structure Fetcher

**Files:**
- Create: `scripts/protein_structures.py`
- Create: `tests/test_protein_structures.py`
- Create: `tests/fixtures/sample_pdb_search.json`
- Create: `tests/fixtures/sample_alphafold.json`

- [ ] **Step 1: Create test fixtures**

Write `tests/fixtures/sample_pdb_search.json` — a realistic PDB search API response for EGFR:

```json
{
  "result_set": [
    {
      "identifier": "1M17",
      "score": 1.0,
      "services": [
        {
          "nodes": [
            {
              "match_context": [{"sequence_identity": 1.0}]
            }
          ]
        }
      ]
    }
  ]
}
```

Write `tests/fixtures/sample_alphafold.json` — an AlphaFold API response:

```json
[
  {
    "entryId": "AF-P00533-F1",
    "gene": "EGFR",
    "uniprotAccession": "P00533",
    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P00533-F1-model_v4.pdb",
    "paeImageUrl": "https://alphafold.ebi.ac.uk/files/AF-P00533-F1-predicted_aligned_error_v4.png"
  }
]
```

- [ ] **Step 2: Write failing tests**

Write `tests/test_protein_structures.py`:

```python
import json
import os
import pytest
from unittest.mock import patch, MagicMock
from scripts.db import init_db
from scripts.protein_structures import (
    search_pdb_for_gene,
    parse_pdb_search_response,
    fetch_alphafold_entry,
    parse_alphafold_response,
    get_uniprot_id,
    DOCKABLE_TARGETS,
    fetch_structure_for_gene,
    store_protein_structure,
)


FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


class TestDockableTargets:
    def test_dockable_targets_is_subset_of_curated(self):
        """DOCKABLE_TARGETS should only contain genes that are realistic docking targets."""
        assert "EGFR" in DOCKABLE_TARGETS
        assert "MTOR" in DOCKABLE_TARGETS
        # Structural proteins should NOT be dockable
        assert "COL1A1" not in DOCKABLE_TARGETS
        assert "COL1A2" not in DOCKABLE_TARGETS

    def test_dockable_targets_has_uniprot_ids(self):
        """Each dockable target should have a UniProt ID for AlphaFold fallback."""
        for gene, info in DOCKABLE_TARGETS.items():
            assert "uniprot" in info, f"{gene} missing uniprot ID"


class TestPDBSearch:
    def test_parse_pdb_search_response(self):
        with open(os.path.join(FIXTURES_DIR, "sample_pdb_search.json")) as f:
            data = json.load(f)
        pdb_ids = parse_pdb_search_response(data)
        assert "1M17" in pdb_ids

    def test_parse_empty_response(self):
        pdb_ids = parse_pdb_search_response({"result_set": []})
        assert pdb_ids == []

    def test_parse_none_response(self):
        pdb_ids = parse_pdb_search_response(None)
        assert pdb_ids == []


class TestAlphaFold:
    def test_parse_alphafold_response(self):
        with open(os.path.join(FIXTURES_DIR, "sample_alphafold.json")) as f:
            data = json.load(f)
        entry = parse_alphafold_response(data)
        assert entry["entry_id"] == "AF-P00533-F1"
        assert "pdb_url" in entry

    def test_parse_empty_alphafold(self):
        entry = parse_alphafold_response([])
        assert entry is None


class TestStoreProteinStructure:
    def test_insert_and_retrieve(self, db_conn):
        init_db(db_conn)
        store_protein_structure(db_conn, {
            "gene_symbol": "EGFR",
            "pdb_id": "1M17",
            "alphafold_id": None,
            "structure_source": "pdb_experimental",
            "resolution_angstroms": 2.6,
            "file_path": "data/structures/proteins/EGFR_1M17.pdb",
        })
        row = db_conn.execute(
            "SELECT * FROM protein_structures WHERE gene_symbol='EGFR'"
        ).fetchone()
        assert row["pdb_id"] == "1M17"
        assert row["structure_source"] == "pdb_experimental"
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_protein_structures.py -v`
Expected: FAIL — `protein_structures` module doesn't exist.

- [ ] **Step 4: Implement protein_structures.py**

Create `scripts/protein_structures.py`:

```python
"""Fetch 3D protein structures from PDB and AlphaFold for keloid targets."""

import json
import logging
import os
import requests
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

PDB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api"

STRUCTURES_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "structures", "proteins")
RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw")

# Curated targets with drugable binding pockets.
# gene_symbol -> {uniprot: UniProt ID, preferred_pdb: known good PDB ID or None}
# Structural proteins (COL1A1/A2, COL3A1) and secreted ligands (TGFB1/B2, WNT3A,
# VEGFA, IL6) are excluded — they lack conventional small-molecule binding pockets.
DOCKABLE_TARGETS = {
    "EGFR":   {"uniprot": "P00533", "preferred_pdb": "1M17"},
    "MTOR":   {"uniprot": "P42345", "preferred_pdb": "4DRH"},
    "PIK3CA": {"uniprot": "P42336", "preferred_pdb": "4OVU"},
    "AKT1":   {"uniprot": "P31749", "preferred_pdb": "3O96"},
    "JAK2":   {"uniprot": "O60674", "preferred_pdb": "3FUP"},
    "STAT3":  {"uniprot": "P40763", "preferred_pdb": "6NUQ"},
    "MAPK1":  {"uniprot": "P28482", "preferred_pdb": "4QTA"},
    "MAPK3":  {"uniprot": "P27361", "preferred_pdb": "4QTB"},
    "KDR":    {"uniprot": "P35968", "preferred_pdb": "1YWN"},
    "CYP24A1":{"uniprot": "Q07973", "preferred_pdb": "3K9V"},
    "SMAD3":  {"uniprot": "P84022", "preferred_pdb": None},
    "TNF":    {"uniprot": "P01375", "preferred_pdb": "2AZ5"},
}


def search_pdb_for_gene(gene_symbol, uniprot_id):
    """Search RCSB PDB for crystal structures of a protein by UniProt ID.
    Returns list of PDB IDs sorted by resolution (best first)."""
    cache_path = os.path.join(RAW_DIR, f"pdb_search_{gene_symbol}.json")
    if os.path.exists(cache_path):
        logger.info(f"Using cached PDB search for {gene_symbol}")
        with open(cache_path) as f:
            return parse_pdb_search_response(json.load(f))

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            "paginate": {"start": 0, "rows": 5},
        },
    }

    try:
        resp = requests.post(PDB_SEARCH_URL, json=query, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except (requests.RequestException, ValueError) as e:
        logger.warning(f"PDB search failed for {gene_symbol}: {e}")
        return []

    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)

    return parse_pdb_search_response(data)


def parse_pdb_search_response(data):
    """Extract PDB IDs from RCSB search API response."""
    if not data or "result_set" not in data:
        return []
    return [entry["identifier"] for entry in data.get("result_set", [])]


def download_pdb_file(pdb_id, gene_symbol):
    """Download a .pdb file from RCSB. Returns local file path."""
    os.makedirs(STRUCTURES_DIR, exist_ok=True)
    file_path = os.path.join(STRUCTURES_DIR, f"{gene_symbol}_{pdb_id}.pdb")
    if os.path.exists(file_path):
        logger.info(f"PDB file already downloaded: {file_path}")
        return file_path

    url = f"{PDB_DOWNLOAD_URL}/{pdb_id}.pdb"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    with open(file_path, "w") as f:
        f.write(resp.text)
    logger.info(f"Downloaded {pdb_id}.pdb for {gene_symbol}")
    return file_path


def fetch_alphafold_entry(uniprot_id):
    """Fetch AlphaFold prediction entry for a UniProt ID."""
    cache_path = os.path.join(RAW_DIR, f"alphafold_{uniprot_id}.json")
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            return parse_alphafold_response(json.load(f))

    url = f"{ALPHAFOLD_API_URL}/prediction/{uniprot_id}"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except (requests.RequestException, ValueError) as e:
        logger.warning(f"AlphaFold lookup failed for {uniprot_id}: {e}")
        return None

    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)

    return parse_alphafold_response(data)


def parse_alphafold_response(data):
    """Parse AlphaFold API response into a dict with entry_id and pdb_url."""
    if not data:
        return None
    # AlphaFold API returns a list (one entry per UniProt ID)
    entry = data[0] if isinstance(data, list) else data
    return {
        "entry_id": entry.get("entryId"),
        "pdb_url": entry.get("pdbUrl"),
        "uniprot": entry.get("uniprotAccession"),
    }


def download_alphafold_pdb(entry, gene_symbol):
    """Download AlphaFold predicted structure. Returns local file path."""
    os.makedirs(STRUCTURES_DIR, exist_ok=True)
    af_id = entry["entry_id"]
    file_path = os.path.join(STRUCTURES_DIR, f"{gene_symbol}_{af_id}.pdb")
    if os.path.exists(file_path):
        return file_path

    resp = requests.get(entry["pdb_url"], timeout=60)
    resp.raise_for_status()

    with open(file_path, "w") as f:
        f.write(resp.text)
    logger.info(f"Downloaded AlphaFold structure for {gene_symbol}")
    return file_path


def get_uniprot_id(gene_symbol):
    """Look up UniProt ID for a gene from DOCKABLE_TARGETS."""
    info = DOCKABLE_TARGETS.get(gene_symbol)
    return info["uniprot"] if info else None


def fetch_structure_for_gene(gene_symbol):
    """Fetch the best available 3D structure for a gene.
    Priority: preferred PDB > PDB search > AlphaFold.
    Returns dict with structure metadata or None."""
    info = DOCKABLE_TARGETS.get(gene_symbol)
    if not info:
        logger.info(f"Skipping {gene_symbol} — not in DOCKABLE_TARGETS")
        return None

    uniprot_id = info["uniprot"]

    # Try preferred PDB first
    if info.get("preferred_pdb"):
        pdb_id = info["preferred_pdb"]
        try:
            file_path = download_pdb_file(pdb_id, gene_symbol)
            return {
                "gene_symbol": gene_symbol,
                "pdb_id": pdb_id,
                "alphafold_id": None,
                "structure_source": "pdb_experimental",
                "resolution_angstroms": None,
                "file_path": file_path,
            }
        except requests.RequestException:
            logger.warning(f"Failed to download preferred PDB {pdb_id} for {gene_symbol}")

    # Try PDB search
    pdb_ids = search_pdb_for_gene(gene_symbol, uniprot_id)
    if pdb_ids:
        pdb_id = pdb_ids[0]
        try:
            file_path = download_pdb_file(pdb_id, gene_symbol)
            return {
                "gene_symbol": gene_symbol,
                "pdb_id": pdb_id,
                "alphafold_id": None,
                "structure_source": "pdb_experimental",
                "resolution_angstroms": None,
                "file_path": file_path,
            }
        except requests.RequestException:
            logger.warning(f"Failed to download PDB {pdb_id} for {gene_symbol}")

    # Fallback to AlphaFold
    af_entry = fetch_alphafold_entry(uniprot_id)
    if af_entry:
        try:
            file_path = download_alphafold_pdb(af_entry, gene_symbol)
            return {
                "gene_symbol": gene_symbol,
                "pdb_id": None,
                "alphafold_id": af_entry["entry_id"],
                "structure_source": "alphafold_predicted",
                "resolution_angstroms": None,
                "file_path": file_path,
            }
        except requests.RequestException:
            logger.warning(f"Failed to download AlphaFold structure for {gene_symbol}")

    logger.warning(f"No structure found for {gene_symbol}")
    return None


def store_protein_structure(conn, structure_data):
    """Insert a protein structure record into the database."""
    conn.execute(
        """INSERT OR REPLACE INTO protein_structures
           (gene_symbol, pdb_id, alphafold_id, structure_source,
            resolution_angstroms, file_path)
           VALUES (?, ?, ?, ?, ?, ?)""",
        (
            structure_data["gene_symbol"],
            structure_data.get("pdb_id"),
            structure_data.get("alphafold_id"),
            structure_data["structure_source"],
            structure_data.get("resolution_angstroms"),
            structure_data["file_path"],
        ),
    )
    conn.commit()


def main():
    """Fetch structures for all dockable targets and store in DB."""
    conn = get_connection()
    init_db(conn)

    fetched = 0
    for gene_symbol in DOCKABLE_TARGETS:
        structure = fetch_structure_for_gene(gene_symbol)
        if structure:
            store_protein_structure(conn, structure)
            fetched += 1
            logger.info(f"Stored structure for {gene_symbol}: {structure['structure_source']}")

    logger.info(f"Fetched {fetched}/{len(DOCKABLE_TARGETS)} protein structures")
    conn.close()
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_protein_structures.py -v`
Expected: All PASS.

- [ ] **Step 6: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: All PASS.

- [ ] **Step 7: Commit**

```bash
git add scripts/protein_structures.py tests/test_protein_structures.py tests/fixtures/sample_pdb_search.json tests/fixtures/sample_alphafold.json
git commit -m "feat: add protein structure fetcher (PDB + AlphaFold fallback)"
```

---

### Task 4: Drug Structure Fetcher

**Files:**
- Create: `scripts/drug_structures.py`
- Create: `tests/test_drug_structures.py`
- Create: `tests/fixtures/sample_pubchem_3d.json`

- [ ] **Step 1: Create test fixture**

Write `tests/fixtures/sample_pubchem_3d.json`:

```json
{
  "PC_Compounds": [
    {
      "id": {"id": {"cid": 2244}},
      "props": [
        {
          "urn": {"label": "SMILES", "name": "Canonical"},
          "value": {"sval": "CC(=O)OC1=CC=CC=C1C(=O)O"}
        }
      ]
    }
  ]
}
```

- [ ] **Step 2: Write failing tests**

Write `tests/test_drug_structures.py`:

```python
import json
import os
import pytest
from unittest.mock import patch, MagicMock
from scripts.db import init_db
from scripts.drug_structures import (
    search_pubchem_by_name,
    parse_pubchem_cid_response,
    fetch_3d_sdf,
    store_drug_structure,
    convert_sdf_to_pdbqt,
)


FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


class TestPubChemSearch:
    def test_parse_cid_response(self):
        data = {"IdentifierList": {"CID": [2244]}}
        cid = parse_pubchem_cid_response(data)
        assert cid == 2244

    def test_parse_empty_response(self):
        cid = parse_pubchem_cid_response({"IdentifierList": {"CID": []}})
        assert cid is None

    def test_parse_none_response(self):
        cid = parse_pubchem_cid_response(None)
        assert cid is None


class TestStoreDrugStructure:
    def test_insert_and_retrieve(self, db_conn):
        init_db(db_conn)
        store_drug_structure(db_conn, {
            "drug_name": "ASPIRIN",
            "pubchem_cid": "2244",
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "sdf_path": "data/structures/ligands/aspirin.sdf",
            "pdbqt_path": "data/structures/ligands/aspirin.pdbqt",
        })
        row = db_conn.execute(
            "SELECT * FROM drug_structures WHERE drug_name='ASPIRIN'"
        ).fetchone()
        assert row["pubchem_cid"] == "2244"
        assert row["smiles"] == "CC(=O)OC1=CC=CC=C1C(=O)O"


class TestConvertSdfToPdbqt:
    def test_returns_none_when_obabel_missing(self, tmp_path):
        sdf_path = tmp_path / "test.sdf"
        sdf_path.write_text("fake sdf content")
        with patch("scripts.drug_structures.shutil.which", return_value=None):
            result = convert_sdf_to_pdbqt(str(sdf_path))
            assert result is None
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_drug_structures.py -v`
Expected: FAIL — module doesn't exist.

- [ ] **Step 4: Implement drug_structures.py**

Create `scripts/drug_structures.py`:

```python
"""Fetch 3D drug structures from PubChem and convert to PDBQT for docking."""

import json
import logging
import os
import shutil
import subprocess
import requests
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
LIGANDS_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "structures", "ligands")
RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw")


def search_pubchem_by_name(drug_name):
    """Search PubChem for a drug by name. Returns CID or None."""
    cache_path = os.path.join(RAW_DIR, f"pubchem_cid_{drug_name.replace(' ', '_')}.json")
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            return parse_pubchem_cid_response(json.load(f))

    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(drug_name)}/cids/JSON"
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        data = resp.json()
    except (requests.RequestException, ValueError) as e:
        logger.warning(f"PubChem CID lookup failed for {drug_name}: {e}")
        return None

    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)

    return parse_pubchem_cid_response(data)


def parse_pubchem_cid_response(data):
    """Extract first CID from PubChem name search response."""
    if not data:
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


def fetch_3d_sdf(cid, drug_name):
    """Download 3D SDF file from PubChem for a given CID. Returns file path."""
    os.makedirs(LIGANDS_DIR, exist_ok=True)
    safe_name = drug_name.replace(" ", "_").replace("/", "_")
    sdf_path = os.path.join(LIGANDS_DIR, f"{safe_name}.sdf")
    if os.path.exists(sdf_path):
        logger.info(f"SDF already exists for {drug_name}")
        return sdf_path

    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type=3d"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
    except requests.RequestException as e:
        logger.warning(f"Failed to download 3D SDF for {drug_name} (CID {cid}): {e}")
        return None

    with open(sdf_path, "w") as f:
        f.write(resp.text)
    logger.info(f"Downloaded 3D SDF for {drug_name}")
    return sdf_path


def get_smiles(cid):
    """Fetch canonical SMILES from PubChem for a CID."""
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        return props[0]["CanonicalSMILES"] if props else None
    except (requests.RequestException, ValueError, KeyError):
        return None


def convert_sdf_to_pdbqt(sdf_path):
    """Convert SDF to PDBQT format using Open Babel. Returns PDBQT path or None."""
    if not shutil.which("obabel"):
        logger.warning("obabel not found — skipping PDBQT conversion")
        return None

    pdbqt_path = sdf_path.replace(".sdf", ".pdbqt")
    if os.path.exists(pdbqt_path):
        return pdbqt_path

    try:
        subprocess.run(
            ["obabel", sdf_path, "-O", pdbqt_path, "--gen3d", "-h"],
            capture_output=True, text=True, check=True, timeout=60,
        )
        logger.info(f"Converted {os.path.basename(sdf_path)} to PDBQT")
        return pdbqt_path
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        logger.warning(f"PDBQT conversion failed for {sdf_path}: {e}")
        return None


def store_drug_structure(conn, structure_data):
    """Insert a drug structure record into the database."""
    conn.execute(
        """INSERT OR REPLACE INTO drug_structures
           (drug_name, pubchem_cid, smiles, sdf_path, pdbqt_path)
           VALUES (?, ?, ?, ?, ?)""",
        (
            structure_data["drug_name"],
            structure_data.get("pubchem_cid"),
            structure_data.get("smiles"),
            structure_data.get("sdf_path"),
            structure_data.get("pdbqt_path"),
        ),
    )
    conn.commit()


def fetch_structure_for_drug(drug_name):
    """Fetch 3D structure for a drug from PubChem, convert to PDBQT.
    Returns dict with structure metadata or None."""
    cid = search_pubchem_by_name(drug_name)
    if not cid:
        return None

    sdf_path = fetch_3d_sdf(cid, drug_name)
    if not sdf_path:
        return None

    smiles = get_smiles(cid)
    pdbqt_path = convert_sdf_to_pdbqt(sdf_path)

    return {
        "drug_name": drug_name,
        "pubchem_cid": str(cid),
        "smiles": smiles,
        "sdf_path": sdf_path,
        "pdbqt_path": pdbqt_path,
    }


def main(drug_names=None):
    """Fetch structures for specified drugs (or top green/yellow tier drugs)."""
    conn = get_connection()
    init_db(conn)

    if drug_names is None:
        # Load from scores — green and yellow tier drugs
        import json
        scores_path = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "scores.json")
        with open(scores_path) as f:
            scores = json.load(f)
        drug_names = [
            s["drug_name"] for s in scores
            if s.get("severity_tier") in ("green", "yellow")
        ]

    fetched = 0
    for name in drug_names:
        structure = fetch_structure_for_drug(name)
        if structure:
            store_drug_structure(conn, structure)
            fetched += 1

    logger.info(f"Fetched {fetched}/{len(drug_names)} drug structures")
    conn.close()
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_drug_structures.py -v`
Expected: All PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/drug_structures.py tests/test_drug_structures.py tests/fixtures/sample_pubchem_3d.json
git commit -m "feat: add drug structure fetcher (PubChem + PDBQT conversion)"
```

---

### Task 5: Protein Preparation for Docking

**Files:**
- Create: `scripts/protein_prep.py`
- Create: `tests/test_protein_prep.py`

Proteins from PDB contain water molecules, alternate conformations, and lack hydrogens. They need cleaning before Vina can use them.

- [ ] **Step 1: Write failing tests**

Write `tests/test_protein_prep.py`:

```python
import os
import pytest
from scripts.protein_prep import (
    clean_pdb_for_docking,
    convert_pdb_to_pdbqt,
    extract_binding_site_from_hetatm,
)


# Minimal PDB content for testing
SAMPLE_PDB = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      27.340  24.430   2.614  1.00  9.67           N
ATOM      2  CA  ALA A   1      26.266  25.413   2.842  1.00 10.38           C
HETATM    3  C1  ERL A 999      30.000  25.000   3.000  1.00  5.00           C
HETATM    4  O   HOH A 500      10.000  10.000  10.000  1.00 20.00           O
END
"""


class TestCleanPdb:
    def test_removes_water_molecules(self, tmp_path):
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text(SAMPLE_PDB)
        cleaned = clean_pdb_for_docking(str(pdb_file))
        with open(cleaned) as f:
            content = f.read()
        assert "HOH" not in content

    def test_keeps_protein_atoms(self, tmp_path):
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text(SAMPLE_PDB)
        cleaned = clean_pdb_for_docking(str(pdb_file))
        with open(cleaned) as f:
            content = f.read()
        assert "ALA" in content


class TestExtractBindingSite:
    def test_extracts_hetatm_center(self, tmp_path):
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text(SAMPLE_PDB)
        center = extract_binding_site_from_hetatm(str(pdb_file))
        # Should find the non-water HETATM (ERL) at ~(30, 25, 3)
        assert center is not None
        assert abs(center[0] - 30.0) < 1.0
        assert abs(center[1] - 25.0) < 1.0
        assert abs(center[2] - 3.0) < 1.0

    def test_returns_none_when_no_hetatm(self, tmp_path):
        pdb_file = tmp_path / "noligand.pdb"
        pdb_file.write_text("ATOM      1  N   ALA A   1      27.340  24.430   2.614  1.00  9.67           N\nEND\n")
        center = extract_binding_site_from_hetatm(str(pdb_file))
        assert center is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_protein_prep.py -v`
Expected: FAIL — module doesn't exist.

- [ ] **Step 3: Implement protein_prep.py**

Create `scripts/protein_prep.py`:

```python
"""Prepare protein PDB files for docking: clean, add hydrogens, find binding site."""

import logging
import os
import shutil
import subprocess

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Residue names to strip (water, common buffer molecules)
STRIP_RESIDUES = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "ACT", "DMS"}


def clean_pdb_for_docking(pdb_path):
    """Remove water molecules and common crystallization artifacts.
    Writes cleaned file to same directory with _clean suffix.
    Returns path to cleaned file."""
    clean_path = pdb_path.replace(".pdb", "_clean.pdb")
    if os.path.exists(clean_path):
        return clean_path

    with open(pdb_path) as f:
        lines = f.readlines()

    cleaned = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname in STRIP_RESIDUES:
                continue
        cleaned.append(line)

    with open(clean_path, "w") as f:
        f.writelines(cleaned)

    logger.info(f"Cleaned PDB: {os.path.basename(clean_path)}")
    return clean_path


def convert_pdb_to_pdbqt(pdb_path):
    """Convert cleaned PDB to PDBQT format using Open Babel.
    Adds hydrogens and computes partial charges.
    Returns PDBQT path or None if obabel unavailable."""
    if not shutil.which("obabel"):
        logger.warning("obabel not found — cannot convert protein to PDBQT")
        return None

    pdbqt_path = pdb_path.replace(".pdb", ".pdbqt")
    if os.path.exists(pdbqt_path):
        return pdbqt_path

    try:
        subprocess.run(
            ["obabel", pdb_path, "-O", pdbqt_path, "-xr", "-h"],
            capture_output=True, text=True, check=True, timeout=120,
        )
        logger.info(f"Converted to PDBQT: {os.path.basename(pdbqt_path)}")
        return pdbqt_path
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        logger.warning(f"Protein PDBQT conversion failed: {e}")
        return None


def extract_binding_site_from_hetatm(pdb_path):
    """Find the centroid of non-water HETATM records as a proxy for the binding site.
    Returns (x, y, z) tuple or None if no ligand found."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname in STRIP_RESIDUES:
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    continue

    if not coords:
        return None

    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    return (round(cx, 3), round(cy, 3), round(cz, 3))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_protein_prep.py -v`
Expected: All PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/protein_prep.py tests/test_protein_prep.py
git commit -m "feat: add protein preparation module (clean PDB, find binding site)"
```

---

### Task 6: AutoDock Vina Docking Module

**Files:**
- Create: `scripts/docking.py`
- Create: `tests/test_docking.py`
- Create: `tests/fixtures/sample_vina_output.txt`

- [ ] **Step 1: Create Vina output fixture**

Write `tests/fixtures/sample_vina_output.txt`:

```
AutoDock Vina v1.2.5
Scoring function : vina

Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 42
Performing search ... done.
Refining results ... done.

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -8.300      0.000      0.000
   2       -7.900      1.234      2.456
   3       -7.200      2.345      3.567
   4       -6.800      3.456      4.678
   5       -6.500      4.567      5.789
```

- [ ] **Step 2: Write failing tests**

Write `tests/test_docking.py`:

```python
import os
import pytest
from unittest.mock import patch, MagicMock
from scripts.db import init_db
from scripts.docking import (
    parse_vina_output,
    build_vina_command,
    store_docking_result,
    get_docking_pairs,
)


FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


class TestParseVinaOutput:
    def test_parses_binding_energies(self):
        with open(os.path.join(FIXTURES_DIR, "sample_vina_output.txt")) as f:
            output = f.read()
        results = parse_vina_output(output)
        assert results["best_energy"] == -8.3
        assert results["num_poses"] == 5

    def test_parses_empty_output(self):
        results = parse_vina_output("")
        assert results["best_energy"] is None
        assert results["num_poses"] == 0


class TestBuildVinaCommand:
    def test_builds_correct_command(self):
        cmd = build_vina_command(
            receptor="protein.pdbqt",
            ligand="drug.pdbqt",
            center=(30.0, 25.0, 3.0),
            box_size=22.0,
            output="output.pdbqt",
            exhaustiveness=8,
        )
        assert cmd[0] == "vina"
        assert "--receptor" in cmd
        assert "protein.pdbqt" in cmd
        assert "--ligand" in cmd
        assert "drug.pdbqt" in cmd
        assert "--center_x" in cmd
        assert "--exhaustiveness" in cmd


class TestGetDockingPairs:
    def test_returns_pairs_for_matching_drugs_and_proteins(self, db_conn):
        init_db(db_conn)
        # Insert test data
        db_conn.execute(
            """INSERT INTO keloid_targets (gene_symbol, target_name, pathway, evidence_type, evidence_strength, source)
               VALUES ('EGFR', 'EGF receptor', 'EGFR', 'expression', 0.8, 'manual')"""
        )
        db_conn.execute(
            "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (1, 'IMATINIB', 'imatinib', 'approved')"
        )
        db_conn.execute(
            "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (1, 'EGFR', 'inhibitor', 'DGIdb')"
        )
        db_conn.execute(
            """INSERT INTO protein_structures (id, gene_symbol, pdb_id, structure_source, file_path)
               VALUES (1, 'EGFR', '1M17', 'pdb_experimental', 'proteins/EGFR.pdbqt')"""
        )
        db_conn.execute(
            """INSERT INTO drug_structures (id, drug_name, pubchem_cid, smiles, sdf_path, pdbqt_path)
               VALUES (1, 'IMATINIB', '5291', 'C', 'ligands/imatinib.sdf', 'ligands/imatinib.pdbqt')"""
        )
        db_conn.commit()

        pairs = get_docking_pairs(db_conn)
        assert len(pairs) == 1
        assert pairs[0]["drug_name"] == "IMATINIB"
        assert pairs[0]["gene_symbol"] == "EGFR"

    def test_skips_drugs_without_pdbqt(self, db_conn):
        init_db(db_conn)
        db_conn.execute(
            """INSERT INTO keloid_targets (gene_symbol, target_name, pathway, evidence_type, evidence_strength, source)
               VALUES ('EGFR', 'EGF receptor', 'EGFR', 'expression', 0.8, 'manual')"""
        )
        db_conn.execute(
            "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (1, 'NODOCK', 'nodock', 'approved')"
        )
        db_conn.execute(
            "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (1, 'EGFR', 'inhibitor', 'DGIdb')"
        )
        db_conn.execute(
            """INSERT INTO protein_structures (id, gene_symbol, pdb_id, structure_source, file_path)
               VALUES (1, 'EGFR', '1M17', 'pdb_experimental', 'proteins/EGFR.pdbqt')"""
        )
        db_conn.execute(
            """INSERT INTO drug_structures (id, drug_name, pubchem_cid, smiles, sdf_path, pdbqt_path)
               VALUES (1, 'NODOCK', '999', 'C', 'ligands/nodock.sdf', NULL)"""
        )
        db_conn.commit()

        pairs = get_docking_pairs(db_conn)
        assert len(pairs) == 0


class TestStoreDockingResult:
    def test_insert_and_retrieve(self, db_conn):
        init_db(db_conn)
        db_conn.execute(
            """INSERT INTO protein_structures (id, gene_symbol, pdb_id, structure_source, file_path)
               VALUES (1, 'EGFR', '1M17', 'pdb_experimental', 'p.pdb')"""
        )
        db_conn.execute(
            """INSERT INTO drug_structures (id, drug_name, pubchem_cid, smiles, sdf_path, pdbqt_path)
               VALUES (1, 'IMATINIB', '5291', 'C', 'l.sdf', 'l.pdbqt')"""
        )
        db_conn.commit()

        store_docking_result(db_conn, {
            "drug_structure_id": 1,
            "protein_structure_id": 1,
            "method": "vina",
            "binding_energy_kcal": -8.3,
            "num_poses": 5,
            "best_pose_path": "results/IMATINIB_EGFR_pose1.pdbqt",
        })

        row = db_conn.execute("SELECT * FROM docking_results").fetchone()
        assert row["binding_energy_kcal"] == -8.3
        assert row["num_poses"] == 5
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_docking.py -v`
Expected: FAIL — module doesn't exist.

- [ ] **Step 4: Implement docking.py**

Create `scripts/docking.py`:

```python
"""Run AutoDock Vina molecular docking for drug-protein pairs."""

import logging
import os
import re
import shutil
import subprocess
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DOCKING_RESULTS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "data", "structures", "docking_results"
)

DEFAULT_BOX_SIZE = 22.0
DEFAULT_EXHAUSTIVENESS = 8


def build_vina_command(receptor, ligand, center, box_size, output, exhaustiveness=8):
    """Build the vina CLI command as a list of arguments."""
    return [
        "vina",
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(box_size),
        "--size_y", str(box_size),
        "--size_z", str(box_size),
        "--exhaustiveness", str(exhaustiveness),
        "--out", output,
    ]


def parse_vina_output(stdout_text):
    """Parse Vina stdout to extract binding energies.
    Returns dict with best_energy and num_poses."""
    energies = []
    for line in stdout_text.split("\n"):
        # Match lines like "   1       -8.300      0.000      0.000"
        match = re.match(r"\s+\d+\s+([-\d.]+)\s+", line)
        if match:
            energies.append(float(match.group(1)))

    return {
        "best_energy": energies[0] if energies else None,
        "num_poses": len(energies),
    }


def run_vina_docking(receptor_pdbqt, ligand_pdbqt, center, drug_name, gene_symbol,
                     box_size=DEFAULT_BOX_SIZE, exhaustiveness=DEFAULT_EXHAUSTIVENESS):
    """Run Vina docking for one drug-protein pair. Returns result dict or None."""
    if not shutil.which("vina"):
        logger.warning("vina binary not found — cannot run docking")
        return None

    os.makedirs(DOCKING_RESULTS_DIR, exist_ok=True)
    output_path = os.path.join(DOCKING_RESULTS_DIR, f"{drug_name}_{gene_symbol}_poses.pdbqt")

    if os.path.exists(output_path):
        logger.info(f"Docking result already exists: {drug_name} × {gene_symbol}")
        # Re-parse if result file exists but we need the energy
        # (won't have stdout, so skip)
        return {"best_energy": None, "num_poses": 0, "best_pose_path": output_path}

    cmd = build_vina_command(
        receptor=receptor_pdbqt,
        ligand=ligand_pdbqt,
        center=center,
        box_size=box_size,
        output=output_path,
        exhaustiveness=exhaustiveness,
    )

    logger.info(f"Docking {drug_name} into {gene_symbol}...")
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            logger.warning(f"Vina failed for {drug_name} × {gene_symbol}: {result.stderr}")
            return None

        parsed = parse_vina_output(result.stdout)
        parsed["best_pose_path"] = output_path
        logger.info(
            f"  {drug_name} × {gene_symbol}: best energy = {parsed['best_energy']} kcal/mol "
            f"({parsed['num_poses']} poses)"
        )
        return parsed

    except subprocess.TimeoutExpired:
        logger.warning(f"Vina timed out for {drug_name} × {gene_symbol}")
        return None


def get_docking_pairs(conn):
    """Find all drug-protein pairs that are ready for docking.
    A pair is ready when both the drug has a PDBQT file and the protein
    has a structure, and they share a gene_symbol via drug_targets."""
    cursor = conn.execute("""
        SELECT DISTINCT
            ds.id as drug_structure_id,
            ds.drug_name,
            ds.pdbqt_path as ligand_pdbqt,
            ps.id as protein_structure_id,
            ps.gene_symbol,
            ps.file_path as protein_pdb,
            ps.binding_site_center_x,
            ps.binding_site_center_y,
            ps.binding_site_center_z,
            ps.binding_site_size
        FROM drug_structures ds
        JOIN drugs d ON ds.drug_name = d.drug_name
        JOIN drug_targets dt ON d.id = dt.drug_id
        JOIN protein_structures ps ON dt.gene_symbol = ps.gene_symbol
        WHERE ds.pdbqt_path IS NOT NULL
    """)
    return [dict(row) for row in cursor.fetchall()]


def store_docking_result(conn, result_data):
    """Store a docking result in the database."""
    conn.execute(
        """INSERT INTO docking_results
           (drug_structure_id, protein_structure_id, method,
            binding_energy_kcal, num_poses, best_pose_path)
           VALUES (?, ?, ?, ?, ?, ?)""",
        (
            result_data["drug_structure_id"],
            result_data["protein_structure_id"],
            result_data.get("method", "vina"),
            result_data.get("binding_energy_kcal"),
            result_data.get("num_poses"),
            result_data.get("best_pose_path"),
        ),
    )
    conn.commit()
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_docking.py -v`
Expected: All PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/docking.py tests/test_docking.py tests/fixtures/sample_vina_output.txt
git commit -m "feat: add AutoDock Vina docking module"
```

---

### Task 7: Pipeline Entry Point (Stage 06)

**Files:**
- Create: `scripts/06_structural_docking.py`
- Modify: `run_pipeline.sh`

- [ ] **Step 1: Create entry point**

Write `scripts/06_structural_docking.py`:

```python
"""Stage 06: Structural docking — fetch structures and run AutoDock Vina."""

from scripts.protein_structures import main as fetch_proteins
from scripts.drug_structures import main as fetch_drugs
from scripts.protein_prep import clean_pdb_for_docking, convert_pdb_to_pdbqt, extract_binding_site_from_hetatm
from scripts.docking import get_docking_pairs, run_vina_docking, store_docking_result
from scripts.db import get_connection, init_db

import logging
import os

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def prepare_proteins(conn):
    """Clean and convert all fetched protein structures. Update binding site coords."""
    cursor = conn.execute("SELECT id, gene_symbol, file_path FROM protein_structures")
    for row in cursor.fetchall():
        pdb_path = row["file_path"]
        if not os.path.exists(pdb_path):
            logger.warning(f"PDB file missing for {row['gene_symbol']}: {pdb_path}")
            continue

        # Clean PDB
        cleaned = clean_pdb_for_docking(pdb_path)

        # Extract binding site from co-crystallized ligand
        center = extract_binding_site_from_hetatm(pdb_path)
        if center:
            conn.execute(
                """UPDATE protein_structures
                   SET binding_site_center_x=?, binding_site_center_y=?, binding_site_center_z=?,
                       file_path=?
                   WHERE id=?""",
                (center[0], center[1], center[2], cleaned, row["id"]),
            )
            logger.info(f"  {row['gene_symbol']}: binding site at ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
        else:
            # Use protein centroid as fallback (less accurate)
            conn.execute(
                "UPDATE protein_structures SET file_path=? WHERE id=?",
                (cleaned, row["id"]),
            )
            logger.warning(f"  {row['gene_symbol']}: no co-crystallized ligand, using protein center")

        # Convert to PDBQT
        convert_pdb_to_pdbqt(cleaned)

    conn.commit()


def run_all_docking(conn):
    """Run Vina docking for all prepared drug-protein pairs."""
    pairs = get_docking_pairs(conn)
    logger.info(f"Found {len(pairs)} drug-protein pairs to dock")

    completed = 0
    for pair in pairs:
        # Need PDBQT for protein
        protein_pdbqt = pair["protein_pdb"].replace(".pdb", ".pdbqt")
        if not os.path.exists(protein_pdbqt):
            logger.warning(f"Protein PDBQT missing for {pair['gene_symbol']}")
            continue

        center_x = pair.get("binding_site_center_x")
        center_y = pair.get("binding_site_center_y")
        center_z = pair.get("binding_site_center_z")
        if center_x is None:
            logger.warning(f"No binding site for {pair['gene_symbol']} — skipping")
            continue

        center = (center_x, center_y, center_z)
        box_size = pair.get("binding_site_size") or 22.0

        result = run_vina_docking(
            receptor_pdbqt=protein_pdbqt,
            ligand_pdbqt=pair["ligand_pdbqt"],
            center=center,
            drug_name=pair["drug_name"],
            gene_symbol=pair["gene_symbol"],
            box_size=box_size,
        )

        if result and result.get("best_energy") is not None:
            store_docking_result(conn, {
                "drug_structure_id": pair["drug_structure_id"],
                "protein_structure_id": pair["protein_structure_id"],
                "method": "vina",
                "binding_energy_kcal": result["best_energy"],
                "num_poses": result["num_poses"],
                "best_pose_path": result["best_pose_path"],
            })
            completed += 1

    logger.info(f"Completed {completed}/{len(pairs)} docking runs")


def main():
    logger.info("=== Phase 1: Fetching protein structures ===")
    fetch_proteins()

    logger.info("=== Phase 2: Fetching drug structures ===")
    fetch_drugs()

    logger.info("=== Phase 3: Preparing proteins ===")
    conn = get_connection()
    init_db(conn)
    prepare_proteins(conn)

    logger.info("=== Phase 4: Running docking ===")
    run_all_docking(conn)

    # Summary
    cursor = conn.execute("""
        SELECT ds.drug_name, ps.gene_symbol, dr.binding_energy_kcal
        FROM docking_results dr
        JOIN drug_structures ds ON dr.drug_structure_id = ds.id
        JOIN protein_structures ps ON dr.protein_structure_id = ps.id
        ORDER BY dr.binding_energy_kcal ASC
    """)
    results = cursor.fetchall()
    if results:
        print(f"\n{'Drug':<30}{'Target':<12}{'Energy (kcal/mol)'}")
        print("-" * 55)
        for r in results[:20]:
            print(f"{r['drug_name']:<30}{r['gene_symbol']:<12}{r['binding_energy_kcal']:.1f}")

    conn.close()


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Update run_pipeline.sh**

Add stage 06 to `run_pipeline.sh` before the report generation step:

```bash
echo "=== Step 6: Structural docking (this may take hours) ==="
python3 -m scripts.06_structural_docking
```

- [ ] **Step 3: Verify it imports correctly**

```bash
python3 -c "from scripts import protein_structures, drug_structures, docking, protein_prep; print('All modules import OK')"
```

- [ ] **Step 4: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: All tests PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/06_structural_docking.py run_pipeline.sh
git commit -m "feat: add pipeline stage 06 — structural docking orchestration"
```

---

### Task 8: Integrate Docking Scores into Composite Score

**Files:**
- Modify: `scripts/score_candidates.py`
- Modify: `tests/test_scoring.py`

- [ ] **Step 1: Write failing tests for structural score integration**

Add to `tests/test_scoring.py`:

```python
class TestStructuralScoreIntegration:
    def test_compute_structural_binding_score_strong_binder(self):
        """Energy < -7.0 should produce a high structural score."""
        from scripts.score_candidates import compute_structural_binding_score
        score = compute_structural_binding_score(-8.5)
        assert score > 0.7

    def test_compute_structural_binding_score_weak_binder(self):
        """Energy > -5.0 should produce a low structural score."""
        from scripts.score_candidates import compute_structural_binding_score
        score = compute_structural_binding_score(-4.0)
        assert score < 0.3

    def test_compute_structural_binding_score_none(self):
        """No docking data should return 0."""
        from scripts.score_candidates import compute_structural_binding_score
        score = compute_structural_binding_score(None)
        assert score == 0.0

    def test_composite_score_with_structural_data(self, populated_db):
        """When structural data exists, composite score should incorporate it."""
        from scripts.score_candidates import query_overlaps, compute_scores, add_structural_scores
        init_db(populated_db)
        # Add structural docking data for DRUG_A targeting EGFR
        populated_db.execute(
            """INSERT INTO protein_structures (id, gene_symbol, pdb_id, structure_source, file_path)
               VALUES (1, 'EGFR', '1M17', 'pdb_experimental', 'p.pdb')"""
        )
        populated_db.execute(
            """INSERT INTO drug_structures (id, drug_name, pubchem_cid, smiles, sdf_path, pdbqt_path)
               VALUES (1, 'DRUG_A', '123', 'C', 'l.sdf', 'l.pdbqt')"""
        )
        populated_db.execute(
            """INSERT INTO docking_results (drug_structure_id, protein_structure_id, method, binding_energy_kcal, num_poses)
               VALUES (1, 1, 'vina', -8.5, 5)"""
        )
        populated_db.commit()

        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        scores = add_structural_scores(populated_db, scores)

        drug_a = [s for s in scores if s["drug_name"] == "DRUG_A"][0]
        assert "structural_binding_score" in drug_a
        assert drug_a["structural_binding_score"] > 0.0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_scoring.py::TestStructuralScoreIntegration -v`
Expected: FAIL — functions don't exist yet.

- [ ] **Step 3: Add structural scoring functions to score_candidates.py**

Add these functions to `scripts/score_candidates.py`:

```python
def compute_structural_binding_score(binding_energy_kcal):
    """Convert Vina binding energy to a 0.0-1.0 score.
    Scale: -10 kcal/mol → 1.0, -5 kcal/mol → 0.0.
    More negative = stronger binding = higher score."""
    if binding_energy_kcal is None:
        return 0.0
    # Clamp to [-10, -3] range then normalize
    clamped = max(-10.0, min(-3.0, binding_energy_kcal))
    # -10 → 1.0, -3 → 0.0
    return round((clamped - (-3.0)) / (-10.0 - (-3.0)), 4)


def add_structural_scores(conn, scores):
    """Look up docking results and add structural_binding_score to each drug score.
    Uses the best (most negative) binding energy across all targets for each drug."""
    cursor = conn.execute("""
        SELECT ds.drug_name, MIN(dr.binding_energy_kcal) as best_energy
        FROM docking_results dr
        JOIN drug_structures ds ON dr.drug_structure_id = ds.id
        GROUP BY ds.drug_name
    """)
    docking_data = {row["drug_name"]: row["best_energy"] for row in cursor.fetchall()}

    for s in scores:
        energy = docking_data.get(s["drug_name"])
        s["best_binding_energy"] = energy
        s["structural_binding_score"] = compute_structural_binding_score(energy)

    return scores
```

Update the weights at the top of `score_candidates.py` — when structural data is available, redistribute weights:

```python
# Updated weights (applied when structural data exists)
W_OVERLAP_STRUCTURAL = 0.30
W_EVIDENCE_STRUCTURAL = 0.30
W_STRUCTURAL = 0.25
W_MULTITARGET_STRUCTURAL = 0.15
```

Update `compute_scores` to accept an optional `include_structural` flag, or keep the existing scoring untouched and add a separate `recompute_with_structural` function:

```python
def recompute_with_structural(scores):
    """Recompute composite scores incorporating structural binding data.
    Only applied to drugs that have docking results."""
    has_structural = [s for s in scores if s.get("structural_binding_score", 0) > 0]
    no_structural = [s for s in scores if s.get("structural_binding_score", 0) == 0]

    # Recompute scores for drugs with structural data
    if has_structural:
        max_overlap = max(s["pathway_overlap_count"] for s in has_structural) or 1
        for s in has_structural:
            norm_overlap = s["pathway_overlap_count"] / max_overlap
            s["composite_score"] = round(
                norm_overlap * W_OVERLAP_STRUCTURAL
                + s["avg_evidence_strength"] * W_EVIDENCE_STRUCTURAL
                + s["structural_binding_score"] * W_STRUCTURAL
                + s["multi_target_bonus"] * W_MULTITARGET_STRUCTURAL,
                4,
            )

    all_scores = has_structural + no_structural
    all_scores.sort(key=lambda s: s["composite_score"], reverse=True)
    return all_scores
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_scoring.py -v`
Expected: All PASS (old tests + new structural tests).

- [ ] **Step 5: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: All PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/score_candidates.py tests/test_scoring.py
git commit -m "feat: integrate docking binding energy into composite scoring"
```

---

### Task 9: Update .gitignore and Documentation

**Files:**
- Modify: `.gitignore`
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update .gitignore**

Add to `.gitignore`:

```
# Structure files (large, downloaded from APIs)
data/structures/
```

- [ ] **Step 2: Update CLAUDE.md with new module info**

Add to `CLAUDE.md`:

```markdown
## Structural Docking (Stage 06)

New modules for protein-drug molecular docking:
- `scripts/protein_structures.py` — fetch PDB/AlphaFold protein structures
- `scripts/drug_structures.py` — fetch PubChem drug structures, convert to PDBQT
- `scripts/protein_prep.py` — clean PDB files, find binding sites
- `scripts/docking.py` — run AutoDock Vina, parse results

New DB tables: `protein_structures`, `drug_structures`, `docking_results`.

Requires: `brew install open-babel autodock-vina` and `pip3 install meeko`.

Docking is CPU-intensive (~2-10 min per pair). Full batch of ~60 pairs takes hours.
```

- [ ] **Step 3: Commit**

```bash
git add .gitignore CLAUDE.md
git commit -m "docs: update gitignore and CLAUDE.md for structural docking"
```

---

## Dependency Summary

| Package | Install | Purpose |
|---------|---------|---------|
| Open Babel | `brew install open-babel` | Format conversion (PDB/SDF → PDBQT) |
| AutoDock Vina | `brew install autodock-vina` | Molecular docking engine |
| meeko | `pip3 install meeko` | Ligand preparation for Vina |

No conda, GPU, or cloud compute needed. Everything runs locally on the M5 MacBook.
