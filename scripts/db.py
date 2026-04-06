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
            id INTEGER PRIMARY KEY,
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
            id INTEGER PRIMARY KEY,
            drug_name TEXT NOT NULL,
            generic_name TEXT NOT NULL UNIQUE,
            approval_status TEXT NOT NULL DEFAULT 'approved',
            original_indication TEXT,
            mechanism_of_action TEXT,
            created_at TEXT DEFAULT (datetime('now'))
        );

        CREATE TABLE IF NOT EXISTS drug_targets (
            id INTEGER PRIMARY KEY,
            drug_id INTEGER NOT NULL REFERENCES drugs(id),
            gene_symbol TEXT NOT NULL,
            action_type TEXT,
            source TEXT NOT NULL,
            created_at TEXT DEFAULT (datetime('now'))
        );

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
    """)
    conn.commit()
