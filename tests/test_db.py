import pytest
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
