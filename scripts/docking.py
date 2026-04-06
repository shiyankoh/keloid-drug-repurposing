"""AutoDock Vina molecular docking for drug-protein pairs.

Runs Vina via subprocess, parses results, and stores them in the DB.
"""

import logging
import os
import re
import shutil
import subprocess

from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
DOCKING_RESULTS_DIR = os.path.join(DATA_DIR, "structures", "docking_results")


# ---------------------------------------------------------------------------
# Command builder
# ---------------------------------------------------------------------------

def build_vina_command(
    receptor: str,
    ligand: str,
    center: tuple,
    box_size: float,
    output: str,
    exhaustiveness: int = 8,
) -> list[str]:
    """Return the CLI args list for an AutoDock Vina run.

    Args:
        receptor: path to receptor .pdbqt
        ligand: path to ligand .pdbqt
        center: (cx, cy, cz) tuple for the docking box center
        box_size: side length (Å) for all three box dimensions
        output: path for the output poses .pdbqt
        exhaustiveness: Vina search exhaustiveness (default 8)
    """
    cx, cy, cz = center
    return [
        "vina",
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", str(box_size),
        "--size_y", str(box_size),
        "--size_z", str(box_size),
        "--exhaustiveness", str(exhaustiveness),
        "--out", output,
    ]


# ---------------------------------------------------------------------------
# Output parser
# ---------------------------------------------------------------------------

def parse_vina_output(stdout_text: str) -> dict:
    """Parse AutoDock Vina stdout and return summary stats.

    Matches result table lines like:
        '   1       -8.300      0.000      0.000'

    Returns:
        {"best_energy": float or None, "num_poses": int}
    """
    if not stdout_text:
        return {"best_energy": None, "num_poses": 0}

    pattern = re.compile(r"\s+\d+\s+([-\d.]+)\s+")
    matches = pattern.findall(stdout_text)

    if not matches:
        return {"best_energy": None, "num_poses": 0}

    return {
        "best_energy": float(matches[0]),
        "num_poses": len(matches),
    }


# ---------------------------------------------------------------------------
# Docking runner
# ---------------------------------------------------------------------------

def run_vina_docking(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    center: tuple,
    drug_name: str,
    gene_symbol: str,
    box_size: float = 22.0,
    exhaustiveness: int = 8,
) -> dict | None:
    """Run AutoDock Vina for a single drug-protein pair.

    Writes output poses to DOCKING_RESULTS_DIR/{drug_name}_{gene_symbol}_poses.pdbqt.

    Returns a dict with parse_vina_output keys plus 'best_pose_path',
    or None if Vina is not installed or the run fails.
    """
    if shutil.which("vina") is None:
        logger.warning("AutoDock Vina not found in PATH — skipping docking")
        return None

    os.makedirs(DOCKING_RESULTS_DIR, exist_ok=True)
    output_path = os.path.join(
        DOCKING_RESULTS_DIR, f"{drug_name}_{gene_symbol}_poses.pdbqt"
    )

    cmd = build_vina_command(
        receptor=receptor_pdbqt,
        ligand=ligand_pdbqt,
        center=center,
        box_size=box_size,
        output=output_path,
        exhaustiveness=exhaustiveness,
    )

    logger.info(f"Running Vina: {drug_name} × {gene_symbol}")
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
        )
        if proc.returncode != 0:
            logger.warning(
                f"Vina exited with code {proc.returncode} for "
                f"{drug_name}/{gene_symbol}: {proc.stderr[:200]}"
            )
            return None
    except subprocess.TimeoutExpired:
        logger.warning(f"Vina timed out for {drug_name}/{gene_symbol}")
        return None
    except OSError as e:
        logger.warning(f"Vina run failed for {drug_name}/{gene_symbol}: {e}")
        return None

    result = parse_vina_output(proc.stdout)
    result["best_pose_path"] = output_path
    logger.info(
        f"Docking done: {drug_name} × {gene_symbol} — "
        f"best energy {result['best_energy']} kcal/mol, "
        f"{result['num_poses']} poses"
    )
    return result


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def get_docking_pairs(conn) -> list[dict]:
    """Return all drug-protein pairs ready for docking.

    Joins drug_structures → drugs → drug_targets → protein_structures.
    Only returns pairs where drug_structures.pdbqt_path IS NOT NULL.
    """
    sql = """
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
    """
    rows = conn.execute(sql).fetchall()
    return [dict(row) for row in rows]


def store_docking_result(conn, result_data: dict) -> None:
    """Insert a docking result into the docking_results table.

    Args:
        result_data: dict with keys:
            drug_structure_id, protein_structure_id, method,
            binding_energy_kcal, num_poses, best_pose_path
    """
    conn.execute(
        """
        INSERT INTO docking_results
            (drug_structure_id, protein_structure_id, method,
             binding_energy_kcal, num_poses, best_pose_path)
        VALUES (?, ?, ?, ?, ?, ?)
        """,
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
    logger.info(
        f"Stored docking result for drug_structure_id={result_data['drug_structure_id']}, "
        f"protein_structure_id={result_data['protein_structure_id']}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    conn = get_connection()
    init_db(conn)

    pairs = get_docking_pairs(conn)
    logger.info(f"Found {len(pairs)} drug-protein pairs ready for docking")

    for pair in pairs:
        center = (
            pair["binding_site_center_x"],
            pair["binding_site_center_y"],
            pair["binding_site_center_z"],
        )
        box_size = pair["binding_site_size"] or 22.0

        result = run_vina_docking(
            receptor_pdbqt=pair["protein_pdb"],
            ligand_pdbqt=pair["ligand_pdbqt"],
            center=center,
            drug_name=pair["drug_name"],
            gene_symbol=pair["gene_symbol"],
            box_size=box_size,
        )

        if result is not None:
            store_docking_result(conn, {
                "drug_structure_id": pair["drug_structure_id"],
                "protein_structure_id": pair["protein_structure_id"],
                "method": "vina",
                "binding_energy_kcal": result.get("best_energy"),
                "num_poses": result.get("num_poses"),
                "best_pose_path": result.get("best_pose_path"),
            })
        else:
            logger.warning(
                f"No docking result for {pair['drug_name']} × {pair['gene_symbol']}"
            )

    conn.close()
    logger.info("Docking run complete.")


if __name__ == "__main__":
    main()
