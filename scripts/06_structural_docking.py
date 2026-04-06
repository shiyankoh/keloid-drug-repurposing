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
            conn.execute(
                "UPDATE protein_structures SET file_path=? WHERE id=?",
                (cleaned, row["id"]),
            )
            logger.warning(f"  {row['gene_symbol']}: no co-crystallized ligand found")

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
        protein_pdb = pair["protein_pdb"]
        # Replace .pdb extension with .pdbqt (handles both _clean.pdb and .pdb)
        if protein_pdb.endswith(".pdb"):
            protein_pdbqt = protein_pdb[:-4] + ".pdbqt"
        else:
            protein_pdbqt = protein_pdb + ".pdbqt"
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
