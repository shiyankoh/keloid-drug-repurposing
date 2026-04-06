"""Tests for scripts/docking.py"""

import os
import pytest

from scripts.docking import (
    DOCKING_RESULTS_DIR,
    build_vina_command,
    parse_vina_output,
    get_docking_pairs,
    store_docking_result,
)
from scripts.db import init_db

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


# ---------------------------------------------------------------------------
# parse_vina_output
# ---------------------------------------------------------------------------

def test_parse_vina_output_with_fixture():
    fixture_path = os.path.join(FIXTURES_DIR, "sample_vina_output.txt")
    with open(fixture_path) as f:
        text = f.read()
    result = parse_vina_output(text)
    assert result["best_energy"] == pytest.approx(-8.3)
    assert result["num_poses"] == 5


def test_parse_vina_output_empty_string():
    result = parse_vina_output("")
    assert result["best_energy"] is None
    assert result["num_poses"] == 0


def test_parse_vina_output_single_pose():
    text = """
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -9.100      0.000      0.000
"""
    result = parse_vina_output(text)
    assert result["best_energy"] == pytest.approx(-9.1)
    assert result["num_poses"] == 1


# ---------------------------------------------------------------------------
# build_vina_command
# ---------------------------------------------------------------------------

def test_build_vina_command_contains_expected_args():
    cmd = build_vina_command(
        receptor="protein.pdbqt",
        ligand="drug.pdbqt",
        center=(1.0, 2.0, 3.0),
        box_size=22.0,
        output="out.pdbqt",
        exhaustiveness=8,
    )
    assert cmd[0] == "vina"
    assert "--receptor" in cmd
    assert "protein.pdbqt" in cmd
    assert "--ligand" in cmd
    assert "drug.pdbqt" in cmd
    assert "--center_x" in cmd
    assert "--center_y" in cmd
    assert "--center_z" in cmd
    assert "1.0" in cmd
    assert "2.0" in cmd
    assert "3.0" in cmd
    assert "--size_x" in cmd
    assert "--size_y" in cmd
    assert "--size_z" in cmd
    assert "22.0" in cmd
    assert "--exhaustiveness" in cmd
    assert "8" in cmd
    assert "--out" in cmd
    assert "out.pdbqt" in cmd


def test_build_vina_command_default_exhaustiveness():
    cmd = build_vina_command(
        receptor="r.pdbqt",
        ligand="l.pdbqt",
        center=(0.0, 0.0, 0.0),
        box_size=20.0,
        output="o.pdbqt",
    )
    idx = cmd.index("--exhaustiveness")
    assert cmd[idx + 1] == "8"


def test_build_vina_command_custom_exhaustiveness():
    cmd = build_vina_command(
        receptor="r.pdbqt",
        ligand="l.pdbqt",
        center=(0.0, 0.0, 0.0),
        box_size=20.0,
        output="o.pdbqt",
        exhaustiveness=16,
    )
    idx = cmd.index("--exhaustiveness")
    assert cmd[idx + 1] == "16"


# ---------------------------------------------------------------------------
# DOCKING_RESULTS_DIR
# ---------------------------------------------------------------------------

def test_docking_results_dir_path():
    assert "docking_results" in DOCKING_RESULTS_DIR
    assert "structures" in DOCKING_RESULTS_DIR


# ---------------------------------------------------------------------------
# get_docking_pairs
# ---------------------------------------------------------------------------

def _seed_full_pair(conn, drug_name="Dexamethasone", gene_symbol="EGFR"):
    """Insert a complete drug-protein pair ready for docking."""
    conn.execute(
        "INSERT INTO drugs (drug_name, generic_name) VALUES (?, ?)",
        (drug_name, drug_name.lower()),
    )
    drug_id = conn.execute(
        "SELECT id FROM drugs WHERE drug_name = ?", (drug_name,)
    ).fetchone()["id"]

    conn.execute(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (?, ?, ?, ?)",
        (drug_id, gene_symbol, "inhibitor", "test"),
    )
    conn.execute(
        """INSERT INTO protein_structures
           (gene_symbol, structure_source, file_path,
            binding_site_center_x, binding_site_center_y, binding_site_center_z,
            binding_site_size)
           VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (gene_symbol, "pdb", f"data/structures/proteins/{gene_symbol}.pdb",
         10.0, 20.0, 30.0, 22.0),
    )
    conn.execute(
        "INSERT INTO drug_structures (drug_name, pdbqt_path) VALUES (?, ?)",
        (drug_name, f"data/structures/drugs/{drug_name}.pdbqt"),
    )
    conn.commit()


def test_get_docking_pairs_returns_pairs(db_conn):
    init_db(db_conn)
    _seed_full_pair(db_conn)

    pairs = get_docking_pairs(db_conn)
    assert len(pairs) == 1
    pair = pairs[0]
    assert pair["drug_name"] == "Dexamethasone"
    assert pair["gene_symbol"] == "EGFR"
    assert pair["ligand_pdbqt"] == "data/structures/drugs/Dexamethasone.pdbqt"
    assert pair["binding_site_center_x"] == pytest.approx(10.0)
    assert pair["binding_site_center_y"] == pytest.approx(20.0)
    assert pair["binding_site_center_z"] == pytest.approx(30.0)
    assert pair["binding_site_size"] == pytest.approx(22.0)


def test_get_docking_pairs_skips_drug_without_pdbqt(db_conn):
    init_db(db_conn)

    # Insert drug + target + protein, but drug_structure has no pdbqt_path
    conn = db_conn
    conn.execute(
        "INSERT INTO drugs (drug_name, generic_name) VALUES (?, ?)",
        ("Rapamycin", "rapamycin"),
    )
    drug_id = conn.execute(
        "SELECT id FROM drugs WHERE drug_name = 'Rapamycin'"
    ).fetchone()["id"]
    conn.execute(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (?, ?, ?, ?)",
        (drug_id, "MTOR", "inhibitor", "test"),
    )
    conn.execute(
        """INSERT INTO protein_structures
           (gene_symbol, structure_source, file_path,
            binding_site_center_x, binding_site_center_y, binding_site_center_z)
           VALUES (?, ?, ?, ?, ?, ?)""",
        ("MTOR", "pdb", "data/structures/proteins/MTOR.pdb", 5.0, 6.0, 7.0),
    )
    # drug_structure exists but pdbqt_path is NULL
    conn.execute(
        "INSERT INTO drug_structures (drug_name, smiles) VALUES (?, ?)",
        ("Rapamycin", "C[C@@H]1CC"),
    )
    conn.commit()

    pairs = get_docking_pairs(db_conn)
    assert len(pairs) == 0


def test_get_docking_pairs_returns_multiple_pairs(db_conn):
    init_db(db_conn)
    _seed_full_pair(db_conn, drug_name="DrugA", gene_symbol="EGFR")
    _seed_full_pair(db_conn, drug_name="DrugB", gene_symbol="JAK2")

    # Add protein structure for JAK2 (EGFR already seeded)
    # Note: _seed_full_pair creates its own protein structure per call
    pairs = get_docking_pairs(db_conn)
    assert len(pairs) == 2


# ---------------------------------------------------------------------------
# store_docking_result
# ---------------------------------------------------------------------------

def test_store_docking_result_insert_and_retrieve(db_conn):
    init_db(db_conn)
    _seed_full_pair(db_conn)

    ds_id = db_conn.execute(
        "SELECT id FROM drug_structures WHERE drug_name = 'Dexamethasone'"
    ).fetchone()["id"]
    ps_id = db_conn.execute(
        "SELECT id FROM protein_structures WHERE gene_symbol = 'EGFR'"
    ).fetchone()["id"]

    result_data = {
        "drug_structure_id": ds_id,
        "protein_structure_id": ps_id,
        "method": "vina",
        "binding_energy_kcal": -8.3,
        "num_poses": 5,
        "best_pose_path": "data/structures/docking_results/Dexamethasone_EGFR_poses.pdbqt",
    }
    store_docking_result(db_conn, result_data)

    row = db_conn.execute(
        "SELECT * FROM docking_results WHERE drug_structure_id = ?", (ds_id,)
    ).fetchone()

    assert row is not None
    assert row["binding_energy_kcal"] == pytest.approx(-8.3)
    assert row["num_poses"] == 5
    assert row["method"] == "vina"
    assert "Dexamethasone_EGFR_poses" in row["best_pose_path"]


def test_store_docking_result_null_energy(db_conn):
    """Store a failed docking result with None energy."""
    init_db(db_conn)
    _seed_full_pair(db_conn)

    ds_id = db_conn.execute(
        "SELECT id FROM drug_structures WHERE drug_name = 'Dexamethasone'"
    ).fetchone()["id"]
    ps_id = db_conn.execute(
        "SELECT id FROM protein_structures WHERE gene_symbol = 'EGFR'"
    ).fetchone()["id"]

    result_data = {
        "drug_structure_id": ds_id,
        "protein_structure_id": ps_id,
        "method": "vina",
        "binding_energy_kcal": None,
        "num_poses": 0,
        "best_pose_path": None,
    }
    store_docking_result(db_conn, result_data)

    row = db_conn.execute(
        "SELECT * FROM docking_results WHERE drug_structure_id = ?", (ds_id,)
    ).fetchone()

    assert row is not None
    assert row["binding_energy_kcal"] is None
    assert row["num_poses"] == 0
