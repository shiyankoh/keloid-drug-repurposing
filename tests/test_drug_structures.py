"""Tests for scripts/drug_structures.py"""

import json
import os
from unittest.mock import patch, MagicMock

import pytest

from scripts.drug_structures import (
    parse_pubchem_cid_response,
    store_drug_structure,
    convert_sdf_to_pdbqt,
)
from scripts.db import init_db

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


def _load_fixture(name):
    with open(os.path.join(FIXTURES_DIR, name)) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# parse_pubchem_cid_response
# ---------------------------------------------------------------------------

def test_parse_pubchem_cid_response_valid():
    data = _load_fixture("sample_pubchem_3d.json")
    cid = parse_pubchem_cid_response(data)
    assert cid == 2244


def test_parse_pubchem_cid_response_empty_cid_list():
    data = {"IdentifierList": {"CID": []}}
    assert parse_pubchem_cid_response(data) is None


def test_parse_pubchem_cid_response_none_input():
    assert parse_pubchem_cid_response(None) is None


def test_parse_pubchem_cid_response_missing_identifier_list():
    assert parse_pubchem_cid_response({}) is None


# ---------------------------------------------------------------------------
# convert_sdf_to_pdbqt
# ---------------------------------------------------------------------------

def test_convert_sdf_to_pdbqt_returns_none_when_obabel_missing(tmp_path):
    sdf_file = tmp_path / "test.sdf"
    sdf_file.write_text("fake sdf content")

    with patch("scripts.drug_structures.shutil.which", return_value=None):
        result = convert_sdf_to_pdbqt(str(sdf_file))

    assert result is None


def test_convert_sdf_to_pdbqt_returns_path_when_obabel_present(tmp_path):
    sdf_file = tmp_path / "aspirin.sdf"
    sdf_file.write_text("fake sdf content")
    expected_pdbqt = str(sdf_file).replace(".sdf", ".pdbqt")

    with patch("scripts.drug_structures.shutil.which", return_value="/usr/bin/obabel"), \
         patch("scripts.drug_structures.subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        # Create a dummy pdbqt file so the function can return the path
        open(expected_pdbqt, "w").close()
        result = convert_sdf_to_pdbqt(str(sdf_file))

    assert result == expected_pdbqt


# ---------------------------------------------------------------------------
# store_drug_structure
# ---------------------------------------------------------------------------

def test_store_drug_structure_insert_and_retrieve(db_conn):
    init_db(db_conn)

    structure = {
        "drug_name": "Aspirin",
        "pubchem_cid": "2244",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "sdf_path": "data/structures/ligands/Aspirin.sdf",
        "pdbqt_path": "data/structures/ligands/Aspirin.pdbqt",
    }
    store_drug_structure(db_conn, structure)

    row = db_conn.execute(
        "SELECT * FROM drug_structures WHERE drug_name = 'Aspirin'"
    ).fetchone()

    assert row is not None
    assert row["pubchem_cid"] == "2244"
    assert row["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"
    assert row["sdf_path"] == "data/structures/ligands/Aspirin.sdf"
    assert row["pdbqt_path"] == "data/structures/ligands/Aspirin.pdbqt"


def test_store_drug_structure_replace_on_conflict(db_conn):
    init_db(db_conn)

    first = {
        "drug_name": "Ibuprofen",
        "pubchem_cid": "3672",
        "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "sdf_path": "data/structures/ligands/Ibuprofen.sdf",
        "pdbqt_path": None,
    }
    store_drug_structure(db_conn, first)

    second = {
        "drug_name": "Ibuprofen",
        "pubchem_cid": "3672",
        "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "sdf_path": "data/structures/ligands/Ibuprofen.sdf",
        "pdbqt_path": "data/structures/ligands/Ibuprofen.pdbqt",
    }
    store_drug_structure(db_conn, second)

    rows = db_conn.execute(
        "SELECT * FROM drug_structures WHERE drug_name = 'Ibuprofen'"
    ).fetchall()

    assert len(rows) == 1
    assert rows[0]["pdbqt_path"] == "data/structures/ligands/Ibuprofen.pdbqt"


def test_store_drug_structure_null_fields_allowed(db_conn):
    init_db(db_conn)

    structure = {
        "drug_name": "TestDrug",
        "pubchem_cid": None,
        "smiles": None,
        "sdf_path": None,
        "pdbqt_path": None,
    }
    store_drug_structure(db_conn, structure)

    row = db_conn.execute(
        "SELECT * FROM drug_structures WHERE drug_name = 'TestDrug'"
    ).fetchone()

    assert row is not None
    assert row["pubchem_cid"] is None
    assert row["smiles"] is None
