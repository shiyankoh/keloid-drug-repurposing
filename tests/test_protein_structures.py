"""Tests for scripts/protein_structures.py"""

import json
import os
import pytest

from scripts.protein_structures import (
    DOCKABLE_TARGETS,
    parse_pdb_search_response,
    parse_alphafold_response,
    store_protein_structure,
    get_uniprot_id,
)
from scripts.db import init_db

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


# ---------------------------------------------------------------------------
# DOCKABLE_TARGETS dict
# ---------------------------------------------------------------------------

EXPECTED_DOCKABLE = {"EGFR", "MTOR", "PIK3CA", "AKT1", "JAK2", "STAT3",
                     "MAPK1", "MAPK3", "KDR", "CYP24A1", "SMAD3", "TNF"}

EXCLUDED_NON_DOCKABLE = {"COL1A1", "COL1A2", "COL3A1", "TGFB1", "TGFB2",
                          "WNT3A", "VEGFA", "IL6"}


def test_dockable_targets_contains_expected_genes():
    for gene in EXPECTED_DOCKABLE:
        assert gene in DOCKABLE_TARGETS, f"{gene} missing from DOCKABLE_TARGETS"


def test_dockable_targets_excludes_non_dockable():
    for gene in EXCLUDED_NON_DOCKABLE:
        assert gene not in DOCKABLE_TARGETS, f"{gene} should not be in DOCKABLE_TARGETS"


def test_dockable_targets_all_have_uniprot():
    for gene, info in DOCKABLE_TARGETS.items():
        assert "uniprot" in info, f"{gene} missing 'uniprot' key"
        assert info["uniprot"], f"{gene} has empty uniprot ID"


def test_dockable_targets_all_have_preferred_pdb_key():
    """preferred_pdb can be None, but the key must exist."""
    for gene, info in DOCKABLE_TARGETS.items():
        assert "preferred_pdb" in info, f"{gene} missing 'preferred_pdb' key"


# ---------------------------------------------------------------------------
# get_uniprot_id
# ---------------------------------------------------------------------------

def test_get_uniprot_id_known_gene():
    uid = get_uniprot_id("EGFR")
    assert uid == DOCKABLE_TARGETS["EGFR"]["uniprot"]


def test_get_uniprot_id_unknown_gene():
    assert get_uniprot_id("COL1A1") is None


# ---------------------------------------------------------------------------
# parse_pdb_search_response
# ---------------------------------------------------------------------------

def _load_fixture(name):
    with open(os.path.join(FIXTURES_DIR, name)) as f:
        return json.load(f)


def test_parse_pdb_search_response_valid_fixture():
    data = _load_fixture("sample_pdb_search.json")
    pdb_ids = parse_pdb_search_response(data)
    assert pdb_ids == ["1IVO", "2GS6", "3VJN"]


def test_parse_pdb_search_response_empty_result_set():
    data = {"result_set": [], "total_count": 0}
    assert parse_pdb_search_response(data) == []


def test_parse_pdb_search_response_missing_result_set():
    data = {"total_count": 0}
    assert parse_pdb_search_response(data) == []


def test_parse_pdb_search_response_none_input():
    assert parse_pdb_search_response(None) == []


# ---------------------------------------------------------------------------
# parse_alphafold_response
# ---------------------------------------------------------------------------

def test_parse_alphafold_response_valid_fixture():
    data = _load_fixture("sample_alphafold.json")
    result = parse_alphafold_response(data)
    assert result is not None
    assert result["entry_id"] == "AF-P00533-F1"
    assert result["pdb_url"] == "https://alphafold.ebi.ac.uk/files/AF-P00533-F1-model_v4.pdb"
    assert result["uniprot"] == "P00533"


def test_parse_alphafold_response_empty_list():
    assert parse_alphafold_response([]) is None


def test_parse_alphafold_response_none_input():
    assert parse_alphafold_response(None) is None


# ---------------------------------------------------------------------------
# store_protein_structure
# ---------------------------------------------------------------------------

def test_store_protein_structure_insert_and_retrieve(db_conn):
    init_db(db_conn)

    structure = {
        "gene_symbol": "EGFR",
        "pdb_id": "1IVO",
        "alphafold_id": None,
        "structure_source": "pdb",
        "resolution_angstroms": 2.5,
        "file_path": "data/structures/proteins/EGFR_1IVO.pdb",
    }
    store_protein_structure(db_conn, structure)

    row = db_conn.execute(
        "SELECT * FROM protein_structures WHERE gene_symbol = 'EGFR'"
    ).fetchone()

    assert row is not None
    assert row["pdb_id"] == "1IVO"
    assert row["structure_source"] == "pdb"
    assert row["resolution_angstroms"] == 2.5
    assert row["file_path"] == "data/structures/proteins/EGFR_1IVO.pdb"


def test_store_protein_structure_replace_on_conflict(db_conn):
    init_db(db_conn)

    first = {
        "gene_symbol": "MTOR",
        "pdb_id": "1FAP",
        "alphafold_id": None,
        "structure_source": "pdb",
        "resolution_angstroms": 3.0,
        "file_path": "data/structures/proteins/MTOR_1FAP.pdb",
    }
    store_protein_structure(db_conn, first)

    second = {
        "gene_symbol": "MTOR",
        "pdb_id": "4JSP",
        "alphafold_id": None,
        "structure_source": "pdb",
        "resolution_angstroms": 2.1,
        "file_path": "data/structures/proteins/MTOR_4JSP.pdb",
    }
    store_protein_structure(db_conn, second)

    rows = db_conn.execute(
        "SELECT * FROM protein_structures WHERE gene_symbol = 'MTOR'"
    ).fetchall()

    # INSERT OR REPLACE should result in one row (replaced)
    assert len(rows) == 1
    assert rows[0]["pdb_id"] == "4JSP"


def test_store_protein_structure_alphafold_source(db_conn):
    init_db(db_conn)

    structure = {
        "gene_symbol": "SMAD3",
        "pdb_id": None,
        "alphafold_id": "AF-P84022-F1",
        "structure_source": "alphafold",
        "resolution_angstroms": None,
        "file_path": "data/structures/proteins/SMAD3_AF-P84022-F1.pdb",
    }
    store_protein_structure(db_conn, structure)

    row = db_conn.execute(
        "SELECT * FROM protein_structures WHERE gene_symbol = 'SMAD3'"
    ).fetchone()

    assert row is not None
    assert row["pdb_id"] is None
    assert row["alphafold_id"] == "AF-P84022-F1"
    assert row["structure_source"] == "alphafold"
