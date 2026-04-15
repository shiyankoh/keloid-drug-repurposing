import json
import os
import pytest
import sqlite3
import tempfile


def _write_annotations(entries):
    """Write annotation entries to a temp JSON file, return path."""
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False)
    json.dump(entries, f)
    f.close()
    return f.name


@pytest.fixture
def prepped_db():
    """In-memory DB with drug_targets table and a known drug-gene row."""
    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row
    conn.execute("""
        CREATE TABLE drugs (
            id INTEGER PRIMARY KEY,
            drug_name TEXT,
            generic_name TEXT,
            original_indication TEXT,
            mechanism_of_action TEXT,
            approval_status TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE drug_targets (
            id INTEGER PRIMARY KEY,
            drug_id INTEGER,
            gene_symbol TEXT,
            action_type TEXT,
            action_direction TEXT
        )
    """)
    conn.execute("INSERT INTO drugs (id, drug_name) VALUES (1, 'SIROLIMUS')")
    conn.execute("INSERT INTO drug_targets (drug_id, gene_symbol, action_type) VALUES (1, 'MTOR', 'inhibitor')")
    conn.commit()
    return conn


def test_import_happy_path(prepped_db):
    from scripts.import_directionality import import_annotations
    path = _write_annotations([{
        "drug_name": "SIROLIMUS",
        "gene_symbol": "MTOR",
        "action_direction": "inhibitor",
        "confidence": "high",
        "source": "ChEMBL:CHEMBL1614345",
        "keloid_target_role": "pro_keloid",
        "net_effect": "beneficial",
        "notes": "mTOR inhibitor"
    }])
    count = import_annotations(prepped_db, path)
    os.unlink(path)
    assert count == 1
    row = prepped_db.execute(
        "SELECT action_direction FROM drug_targets WHERE gene_symbol = 'MTOR'"
    ).fetchone()
    assert row["action_direction"] == "inhibitor"


def test_import_unknown_drug_skipped(prepped_db):
    from scripts.import_directionality import import_annotations
    path = _write_annotations([{
        "drug_name": "NONEXISTENT_DRUG",
        "gene_symbol": "MTOR",
        "action_direction": "inhibitor",
        "confidence": "low",
        "source": "manual",
        "keloid_target_role": "pro_keloid",
        "net_effect": "beneficial",
        "notes": ""
    }])
    count = import_annotations(prepped_db, path)
    os.unlink(path)
    assert count == 0


def test_import_missing_file_raises():
    from scripts.import_directionality import import_annotations
    conn = sqlite3.connect(":memory:")
    with pytest.raises(FileNotFoundError):
        import_annotations(conn, "/nonexistent/path.json")


def test_import_malformed_entry_raises(prepped_db):
    from scripts.import_directionality import import_annotations
    path = _write_annotations([{"drug_name": "SIROLIMUS"}])  # missing gene_symbol + action_direction
    with pytest.raises(ValueError, match="missing fields"):
        import_annotations(prepped_db, path)
    os.unlink(path)
