"""Tests for scripts/protein_prep.py"""

import os
import pytest
from unittest.mock import patch

from scripts.protein_prep import (
    STRIP_RESIDUES,
    clean_pdb_for_docking,
    convert_pdb_to_pdbqt,
    extract_binding_site_from_hetatm,
)

# ---------------------------------------------------------------------------
# Sample PDB content for fixtures
# ---------------------------------------------------------------------------

SAMPLE_PDB = """\
HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      27.340  24.430   2.614  1.00  9.67           N
ATOM      2  CA  ALA A   1      26.266  25.413   2.842  1.00 10.38           C
HETATM    3  C1  ERL A 999      30.000  25.000   3.000  1.00  5.00           C
HETATM    4  O   HOH A 500      10.000  10.000  10.000  1.00 20.00           O
END
"""

SAMPLE_PDB_NO_HETATM = """\
HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      27.340  24.430   2.614  1.00  9.67           N
ATOM      2  CA  ALA A   1      26.266  25.413   2.842  1.00 10.38           C
END
"""


# ---------------------------------------------------------------------------
# STRIP_RESIDUES
# ---------------------------------------------------------------------------

def test_strip_residues_contains_expected():
    expected = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "ACT", "DMS"}
    assert STRIP_RESIDUES == expected


# ---------------------------------------------------------------------------
# clean_pdb_for_docking
# ---------------------------------------------------------------------------

def test_clean_pdb_removes_water_keeps_protein_and_ligand(tmp_path):
    pdb_file = tmp_path / "test.pdb"
    pdb_file.write_text(SAMPLE_PDB)

    result_path = clean_pdb_for_docking(str(pdb_file))

    assert result_path.endswith("_clean.pdb")
    assert os.path.exists(result_path)

    content = open(result_path).read()

    # Water (HOH) should be removed
    assert "HOH" not in content

    # Protein atoms (ALA) should remain
    assert "ALA" in content

    # Non-water HETATM (ERL) should remain
    assert "ERL" in content


def test_clean_pdb_skips_if_already_exists(tmp_path):
    pdb_file = tmp_path / "test.pdb"
    pdb_file.write_text(SAMPLE_PDB)

    clean_path = str(tmp_path / "test_clean.pdb")
    with open(clean_path, "w") as f:
        f.write("EXISTING CONTENT\n")

    result = clean_pdb_for_docking(str(pdb_file))

    # Should return the existing file without overwriting it
    assert result == clean_path
    assert open(result).read() == "EXISTING CONTENT\n"


def test_clean_pdb_returns_correct_path(tmp_path):
    pdb_file = tmp_path / "myprotein.pdb"
    pdb_file.write_text(SAMPLE_PDB)

    result = clean_pdb_for_docking(str(pdb_file))

    expected = str(tmp_path / "myprotein_clean.pdb")
    assert result == expected


# ---------------------------------------------------------------------------
# extract_binding_site_from_hetatm
# ---------------------------------------------------------------------------

def test_extract_binding_site_finds_centroid(tmp_path):
    pdb_file = tmp_path / "test.pdb"
    pdb_file.write_text(SAMPLE_PDB)

    result = extract_binding_site_from_hetatm(str(pdb_file))

    assert result is not None
    x, y, z = result
    assert abs(x - 30.0) < 0.01
    assert abs(y - 25.0) < 0.01
    assert abs(z - 3.0) < 0.01


def test_extract_binding_site_returns_none_when_no_hetatm(tmp_path):
    pdb_file = tmp_path / "test_no_hetatm.pdb"
    pdb_file.write_text(SAMPLE_PDB_NO_HETATM)

    result = extract_binding_site_from_hetatm(str(pdb_file))

    assert result is None


def test_extract_binding_site_ignores_water_hetatm(tmp_path):
    # File with only HOH HETATM — should return None
    water_only = """\
ATOM      1  N   ALA A   1      27.340  24.430   2.614  1.00  9.67           N
HETATM    2  O   HOH A 500      10.000  10.000  10.000  1.00 20.00           O
END
"""
    pdb_file = tmp_path / "water_only.pdb"
    pdb_file.write_text(water_only)

    result = extract_binding_site_from_hetatm(str(pdb_file))

    assert result is None


# ---------------------------------------------------------------------------
# convert_pdb_to_pdbqt
# ---------------------------------------------------------------------------

def test_convert_pdb_to_pdbqt_returns_none_when_obabel_missing(tmp_path):
    pdb_file = tmp_path / "test.pdb"
    pdb_file.write_text(SAMPLE_PDB)

    with patch("shutil.which", return_value=None):
        result = convert_pdb_to_pdbqt(str(pdb_file))

    assert result is None
