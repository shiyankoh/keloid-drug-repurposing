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
    assert len(cursor.fetchall()) == 6
    conn.close()


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


def test_protein_structures_has_required_columns(db_conn):
    init_db(db_conn)
    cursor = db_conn.execute("PRAGMA table_info(protein_structures)")
    columns = {row["name"]: row["type"] for row in cursor.fetchall()}
    assert "id" in columns
    assert "gene_symbol" in columns
    assert "pdb_id" in columns
    assert "alphafold_id" in columns
    assert "structure_source" in columns
    assert "resolution_angstroms" in columns
    assert "binding_site_center_x" in columns
    assert "binding_site_center_y" in columns
    assert "binding_site_center_z" in columns
    assert "binding_site_size" in columns
    assert "file_path" in columns
    assert "created_at" in columns


def test_drug_structures_has_required_columns(db_conn):
    init_db(db_conn)
    cursor = db_conn.execute("PRAGMA table_info(drug_structures)")
    columns = {row["name"]: row["type"] for row in cursor.fetchall()}
    assert "id" in columns
    assert "drug_name" in columns
    assert "pubchem_cid" in columns
    assert "smiles" in columns
    assert "sdf_path" in columns
    assert "pdbqt_path" in columns
    assert "created_at" in columns


def test_docking_results_has_required_columns(db_conn):
    init_db(db_conn)
    cursor = db_conn.execute("PRAGMA table_info(docking_results)")
    columns = {row["name"]: row["type"] for row in cursor.fetchall()}
    assert "id" in columns
    assert "drug_structure_id" in columns
    assert "protein_structure_id" in columns
    assert "method" in columns
    assert "binding_energy_kcal" in columns
    assert "num_poses" in columns
    assert "best_pose_path" in columns
    assert "created_at" in columns


def test_insert_and_retrieve_protein_structure(db_conn):
    init_db(db_conn)
    db_conn.execute(
        """INSERT INTO protein_structures
           (gene_symbol, pdb_id, alphafold_id, structure_source, resolution_angstroms,
            binding_site_center_x, binding_site_center_y, binding_site_center_z, file_path)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        ("TP53", "1TUP", "AF-P04637-F1", "pdb", 2.15, 25.5, 30.2, 35.1, "/path/to/tp53.pdb"),
    )
    db_conn.commit()

    row = db_conn.execute(
        "SELECT gene_symbol, pdb_id, structure_source FROM protein_structures WHERE gene_symbol = ?"
        , ("TP53",)
    ).fetchone()
    assert row is not None
    assert row["gene_symbol"] == "TP53"
    assert row["pdb_id"] == "1TUP"
    assert row["structure_source"] == "pdb"


def test_docking_results_with_foreign_keys(db_conn):
    import sqlite3 as _sqlite3

    init_db(db_conn)
    db_conn.execute("PRAGMA foreign_keys = ON")

    # Insert a drug structure
    db_conn.execute(
        """INSERT INTO drug_structures (drug_name, smiles, sdf_path, pdbqt_path)
           VALUES (?, ?, ?, ?)""",
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "/path/to/aspirin.sdf", "/path/to/aspirin.pdbqt"),
    )
    db_conn.commit()

    # Insert a protein structure
    db_conn.execute(
        """INSERT INTO protein_structures
           (gene_symbol, pdb_id, structure_source, file_path)
           VALUES (?, ?, ?, ?)""",
        ("TP53", "1TUP", "pdb", "/path/to/tp53.pdb"),
    )
    db_conn.commit()

    # Insert a docking result
    db_conn.execute(
        """INSERT INTO docking_results
           (drug_structure_id, protein_structure_id, method, binding_energy_kcal, num_poses)
           VALUES (?, ?, ?, ?, ?)""",
        (1, 1, "vina", -8.5, 10),
    )
    db_conn.commit()

    # Retrieve and verify
    row = db_conn.execute(
        "SELECT drug_structure_id, protein_structure_id, binding_energy_kcal FROM docking_results WHERE id = ?"
        , (1,)
    ).fetchone()
    assert row is not None
    assert row["drug_structure_id"] == 1
    assert row["protein_structure_id"] == 1
    assert row["binding_energy_kcal"] == -8.5


def test_docking_results_foreign_key_constraint(db_conn):
    import sqlite3 as _sqlite3

    init_db(db_conn)
    db_conn.execute("PRAGMA foreign_keys = ON")

    # Try to insert a docking result with invalid foreign keys
    with pytest.raises(_sqlite3.IntegrityError):
        db_conn.execute(
            """INSERT INTO docking_results
               (drug_structure_id, protein_structure_id, method, binding_energy_kcal)
               VALUES (?, ?, ?, ?)""",
            (999, 999, "vina", -8.5),
        )
        db_conn.commit()
