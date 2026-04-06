"""Protein preparation module for molecular docking.

Cleans PDB files (removes water and crystallographic artifacts) and extracts
binding site coordinates from co-crystallized ligands.
"""

import logging
import os
import shutil
import subprocess

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Residue names to strip from PDB files before docking
STRIP_RESIDUES = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "ACT", "DMS"}


def clean_pdb_for_docking(pdb_path: str) -> str:
    """Remove water and artifact residues from a PDB file.

    Reads pdb_path, strips lines whose residue name (columns 17-20) is in
    STRIP_RESIDUES, and writes the cleaned file to {stem}_clean.pdb.
    Skips if the output file already exists.

    Returns the path to the cleaned PDB file.
    """
    stem = os.path.splitext(pdb_path)[0]
    clean_path = f"{stem}_clean.pdb"

    if os.path.exists(clean_path):
        logger.info(f"Clean PDB already exists, skipping: {clean_path}")
        return clean_path

    kept_lines = []
    with open(pdb_path) as f:
        for line in f:
            # PDB columns 17-20 (0-indexed) hold the residue name
            if line.startswith(("ATOM  ", "HETATM")):
                res_name = line[17:20].strip()
                if res_name in STRIP_RESIDUES:
                    continue
            kept_lines.append(line)

    with open(clean_path, "w") as f:
        f.writelines(kept_lines)

    logger.info(f"Cleaned PDB written to {clean_path}")
    return clean_path


def convert_pdb_to_pdbqt(pdb_path: str) -> str | None:
    """Convert a PDB file to PDBQT format using Open Babel.

    Runs: obabel pdb_path -O output.pdbqt -xr -h

    Returns the path to the PDBQT file, or None if obabel is not found or
    the conversion fails.
    """
    if shutil.which("obabel") is None:
        logger.warning("obabel not found — cannot convert PDB to PDBQT")
        return None

    stem = os.path.splitext(pdb_path)[0]
    pdbqt_path = f"{stem}.pdbqt"

    try:
        subprocess.run(
            ["obabel", pdb_path, "-O", pdbqt_path, "-xr", "-h"],
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        )
    except subprocess.CalledProcessError as e:
        logger.warning(f"obabel conversion failed: {e.stderr}")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("obabel conversion timed out")
        return None

    logger.info(f"Converted {pdb_path} → {pdbqt_path}")
    return pdbqt_path


def extract_binding_site_from_hetatm(pdb_path: str) -> tuple[float, float, float] | None:
    """Find binding site centroid from co-crystallized ligand atoms.

    Reads all HETATM lines that are NOT in STRIP_RESIDUES, extracts x/y/z
    coordinates (PDB columns 30-38, 38-46, 46-54), and returns the centroid
    as a (x, y, z) tuple.

    Returns None if no non-water HETATM atoms are found.
    """
    coords = []

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            res_name = line[17:20].strip()
            if res_name in STRIP_RESIDUES:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                logger.warning(f"Could not parse coordinates from HETATM line: {line.rstrip()}")
                continue

    if not coords:
        logger.info(f"No non-water HETATM atoms found in {pdb_path}")
        return None

    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)

    logger.info(f"Binding site centroid: ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    return (cx, cy, cz)
