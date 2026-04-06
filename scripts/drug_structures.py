"""Fetch 3D drug structures from PubChem and convert to PDBQT.

Pipeline:
  search_pubchem_by_name → fetch_3d_sdf → get_smiles → convert_sdf_to_pdbqt
"""

import json
import logging
import os
import shutil
import subprocess

import requests

from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
LIGANDS_DIR = os.path.join(DATA_DIR, "structures", "ligands")

PUBCHEM_CID_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
PUBCHEM_SDF_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
PUBCHEM_SMILES_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_name(drug_name: str) -> str:
    """Replace spaces and slashes with underscores for use in filenames."""
    return drug_name.replace(" ", "_").replace("/", "_")


# ---------------------------------------------------------------------------
# PubChem search
# ---------------------------------------------------------------------------

def search_pubchem_by_name(drug_name: str) -> int | None:
    """Search PubChem for a drug by name and return its CID.

    Caches the raw API response to data/raw/pubchem_cid_{safe_name}.json.
    Returns CID (int) or None on failure.
    """
    os.makedirs(RAW_DIR, exist_ok=True)
    safe = _safe_name(drug_name)
    cache_path = os.path.join(RAW_DIR, f"pubchem_cid_{safe}.json")

    if os.path.exists(cache_path):
        logger.info(f"Loading cached PubChem CID for {drug_name}")
        with open(cache_path) as f:
            data = json.load(f)
        return parse_pubchem_cid_response(data)

    url = PUBCHEM_CID_URL.format(drug_name=requests.utils.quote(drug_name))
    logger.info(f"Searching PubChem for {drug_name}...")
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        logger.warning(f"PubChem search failed for {drug_name}: {e}")
        return None

    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Cached PubChem CID response for {drug_name}")

    return parse_pubchem_cid_response(data)


def parse_pubchem_cid_response(data) -> int | None:
    """Extract the first CID from a PubChem CID response dict.

    Expected shape: {"IdentifierList": {"CID": [2244, ...]}}
    Returns None for empty/None/malformed input.
    """
    if not data:
        return None
    identifier_list = data.get("IdentifierList")
    if not identifier_list:
        return None
    cids = identifier_list.get("CID", [])
    if not cids:
        return None
    return int(cids[0])


# ---------------------------------------------------------------------------
# 3D SDF download
# ---------------------------------------------------------------------------

def fetch_3d_sdf(cid: int, drug_name: str) -> str | None:
    """Download 3D SDF for a PubChem CID. Skips if file already exists.

    Returns the local file path, or None on failure.
    """
    os.makedirs(LIGANDS_DIR, exist_ok=True)
    safe = _safe_name(drug_name)
    dest = os.path.join(LIGANDS_DIR, f"{safe}.sdf")

    if os.path.exists(dest):
        logger.info(f"SDF already exists: {dest}")
        return dest

    url = PUBCHEM_SDF_URL.format(cid=cid)
    logger.info(f"Downloading 3D SDF for {drug_name} (CID {cid})...")
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException as e:
        logger.warning(f"Failed to download SDF for {drug_name}: {e}")
        return None

    with open(dest, "w") as f:
        f.write(resp.text)
    logger.info(f"Saved SDF to {dest}")
    return dest


# ---------------------------------------------------------------------------
# SMILES
# ---------------------------------------------------------------------------

def get_smiles(cid: int) -> str | None:
    """Fetch canonical SMILES from PubChem for a given CID.

    Returns the SMILES string or None on failure.
    """
    url = PUBCHEM_SMILES_URL.format(cid=cid)
    logger.info(f"Fetching SMILES for CID {cid}...")
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        logger.warning(f"Failed to fetch SMILES for CID {cid}: {e}")
        return None

    try:
        return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except (KeyError, IndexError, TypeError):
        logger.warning(f"Unexpected SMILES response shape for CID {cid}")
        return None


# ---------------------------------------------------------------------------
# PDBQT conversion
# ---------------------------------------------------------------------------

def convert_sdf_to_pdbqt(sdf_path: str) -> str | None:
    """Convert an SDF file to PDBQT using the obabel CLI.

    Runs: obabel input.sdf -O output.pdbqt --gen3d -h
    Returns the PDBQT path on success, or None if obabel is not found.
    """
    if not shutil.which("obabel"):
        logger.warning("obabel not found — skipping PDBQT conversion")
        return None

    pdbqt_path = sdf_path.replace(".sdf", ".pdbqt")
    cmd = ["obabel", sdf_path, "-O", pdbqt_path, "--gen3d", "-h"]
    logger.info(f"Converting {sdf_path} to PDBQT...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            logger.warning(f"obabel returned non-zero exit code: {result.stderr}")
            return None
    except subprocess.TimeoutExpired:
        logger.warning(f"obabel timed out for {sdf_path}")
        return None
    except OSError as e:
        logger.warning(f"obabel subprocess error: {e}")
        return None

    if not os.path.exists(pdbqt_path):
        logger.warning(f"obabel did not produce expected output: {pdbqt_path}")
        return None

    logger.info(f"Saved PDBQT to {pdbqt_path}")
    return pdbqt_path


# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------

def store_drug_structure(conn, structure_data: dict) -> None:
    """INSERT OR REPLACE into drug_structures table.

    Deletes any existing row for drug_name before inserting to maintain
    one row per drug (the table has no UNIQUE constraint on drug_name).
    """
    drug_name = structure_data["drug_name"]
    conn.execute("DELETE FROM drug_structures WHERE drug_name = ?", (drug_name,))
    conn.execute(
        """
        INSERT INTO drug_structures
            (drug_name, pubchem_cid, smiles, sdf_path, pdbqt_path)
        VALUES (?, ?, ?, ?, ?)
        """,
        (
            drug_name,
            structure_data.get("pubchem_cid"),
            structure_data.get("smiles"),
            structure_data.get("sdf_path"),
            structure_data.get("pdbqt_path"),
        ),
    )
    conn.commit()
    logger.info(f"Stored drug structure for {drug_name}")


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def fetch_structure_for_drug(drug_name: str) -> dict | None:
    """Orchestrate PubChem search → SDF download → SMILES → PDBQT conversion.

    Returns a dict with drug_name, pubchem_cid, smiles, sdf_path, pdbqt_path,
    or None if the CID could not be found.
    """
    cid = search_pubchem_by_name(drug_name)
    if cid is None:
        logger.warning(f"No CID found for {drug_name} — skipping")
        return None

    sdf_path = fetch_3d_sdf(cid, drug_name)
    smiles = get_smiles(cid)

    pdbqt_path = None
    if sdf_path:
        pdbqt_path = convert_sdf_to_pdbqt(sdf_path)

    return {
        "drug_name": drug_name,
        "pubchem_cid": str(cid),
        "smiles": smiles,
        "sdf_path": sdf_path,
        "pdbqt_path": pdbqt_path,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(drug_names=None):
    """Fetch and store structures for a list of drugs.

    If drug_names is None, load from data/raw/scores.json and pick green/yellow
    severity tier drugs.
    """
    conn = get_connection()
    init_db(conn)

    if drug_names is None:
        scores_path = os.path.join(DATA_DIR, "raw", "scores.json")
        drug_names = []
        if os.path.exists(scores_path):
            with open(scores_path) as f:
                scores = json.load(f)
            for entry in scores:
                tier = entry.get("severity_tier", "").lower()
                if tier in ("green", "yellow"):
                    name = entry.get("drug_name") or entry.get("generic_name")
                    if name:
                        drug_names.append(name)
            logger.info(f"Loaded {len(drug_names)} green/yellow drugs from scores.json")
        else:
            logger.warning(f"scores.json not found at {scores_path}")

    for drug_name in drug_names:
        structure = fetch_structure_for_drug(drug_name)
        if structure:
            store_drug_structure(conn, structure)
        else:
            logger.warning(f"No structure stored for {drug_name}")

    conn.close()
    logger.info("Drug structure fetch complete.")


if __name__ == "__main__":
    main()
