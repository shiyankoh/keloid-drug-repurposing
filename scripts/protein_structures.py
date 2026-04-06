"""Fetch 3D protein structures for keloid-relevant dockable targets.

Priority: preferred_pdb → PDB search → AlphaFold fallback.
All API responses are cached to data/raw/ as JSON.
"""

import json
import logging
import os

import requests

from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
STRUCTURES_DIR = os.path.join(DATA_DIR, "structures", "proteins")

PDB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"

# Dockable targets only — kinases and enzymes with well-defined binding pockets.
# Excluded: structural proteins (COL1A1/A2, COL3A1) and secreted ligands
# (TGFB1/B2, WNT3A, VEGFA, IL6) which lack tractable small-molecule pockets.
DOCKABLE_TARGETS = {
    "EGFR": {
        "uniprot": "P00533",
        "preferred_pdb": "1IVO",   # EGFR kinase domain, X-ray 2.6 Å
    },
    "MTOR": {
        "uniprot": "P42345",
        "preferred_pdb": "4JSP",   # mTOR kinase domain, X-ray 3.2 Å
    },
    "PIK3CA": {
        "uniprot": "P42336",
        "preferred_pdb": "2RD0",   # PI3Kα, X-ray 3.0 Å
    },
    "AKT1": {
        "uniprot": "P31749",
        "preferred_pdb": "3O96",   # AKT1 kinase domain, X-ray 2.0 Å
    },
    "JAK2": {
        "uniprot": "O60674",
        "preferred_pdb": "3LPB",   # JAK2 kinase domain, X-ray 2.2 Å
    },
    "STAT3": {
        "uniprot": "P40763",
        "preferred_pdb": "1BG1",   # STAT3 SH2 domain, X-ray 2.25 Å
    },
    "MAPK1": {
        "uniprot": "P28482",
        "preferred_pdb": "4GT3",   # ERK2, X-ray 1.8 Å
    },
    "MAPK3": {
        "uniprot": "P27361",
        "preferred_pdb": "2ZOQ",   # ERK1, X-ray 2.4 Å
    },
    "KDR": {
        "uniprot": "P35968",
        "preferred_pdb": "2OH4",   # VEGFR2 kinase domain, X-ray 2.0 Å
    },
    "CYP24A1": {
        "uniprot": "Q07973",
        "preferred_pdb": "3K9V",   # CYP24A1, X-ray 2.8 Å
    },
    "SMAD3": {
        "uniprot": "P84022",
        "preferred_pdb": None,     # no well-curated structure; use AlphaFold
    },
    "TNF": {
        "uniprot": "P01375",
        "preferred_pdb": "2AZ5",   # TNF trimer, X-ray 2.1 Å
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_uniprot_id(gene_symbol: str):
    """Return UniProt ID for a dockable gene, or None if not in the dict."""
    entry = DOCKABLE_TARGETS.get(gene_symbol)
    if entry is None:
        return None
    return entry["uniprot"]


# ---------------------------------------------------------------------------
# PDB search
# ---------------------------------------------------------------------------

def search_pdb_for_gene(gene_symbol: str, uniprot_id: str) -> list[str]:
    """Search RCSB PDB for X-ray crystal structures by UniProt ID.

    Caches the raw API response to data/raw/pdb_search_{gene}.json.
    Returns a list of PDB IDs (may be empty).
    """
    os.makedirs(RAW_DIR, exist_ok=True)
    cache_path = os.path.join(RAW_DIR, f"pdb_search_{gene_symbol}.json")

    if os.path.exists(cache_path):
        logger.info(f"Loading cached PDB search for {gene_symbol}")
        with open(cache_path) as f:
            data = json.load(f)
        return parse_pdb_search_response(data)

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 20},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    }

    logger.info(f"Searching RCSB PDB for {gene_symbol} ({uniprot_id})...")
    try:
        resp = requests.post(PDB_SEARCH_URL, json=query, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        logger.warning(f"PDB search failed for {gene_symbol}: {e}")
        return []

    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Cached PDB search response for {gene_symbol}")

    return parse_pdb_search_response(data)


def parse_pdb_search_response(data) -> list[str]:
    """Extract PDB IDs from RCSB search API response.

    Returns [] for None, empty, or malformed input.
    """
    if not data:
        return []
    result_set = data.get("result_set", [])
    if not result_set:
        return []
    return [entry["identifier"] for entry in result_set if "identifier" in entry]


# ---------------------------------------------------------------------------
# PDB file download
# ---------------------------------------------------------------------------

def download_pdb_file(pdb_id: str, gene_symbol: str) -> str | None:
    """Download a .pdb file from RCSB. Skips if already on disk.

    Returns the local file path, or None on failure.
    """
    os.makedirs(STRUCTURES_DIR, exist_ok=True)
    dest = os.path.join(STRUCTURES_DIR, f"{gene_symbol}_{pdb_id}.pdb")

    if os.path.exists(dest):
        logger.info(f"PDB file already exists: {dest}")
        return dest

    url = PDB_DOWNLOAD_URL.format(pdb_id=pdb_id)
    logger.info(f"Downloading {pdb_id} for {gene_symbol}...")
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException as e:
        logger.warning(f"Failed to download PDB {pdb_id}: {e}")
        return None

    with open(dest, "w") as f:
        f.write(resp.text)
    logger.info(f"Saved {pdb_id} to {dest}")
    return dest


# ---------------------------------------------------------------------------
# AlphaFold
# ---------------------------------------------------------------------------

def fetch_alphafold_entry(uniprot_id: str):
    """Fetch AlphaFold prediction metadata for a UniProt ID.

    Caches to data/raw/alphafold_{uniprot}.json.
    Returns raw API response (list) or None on failure.
    """
    os.makedirs(RAW_DIR, exist_ok=True)
    cache_path = os.path.join(RAW_DIR, f"alphafold_{uniprot_id}.json")

    if os.path.exists(cache_path):
        logger.info(f"Loading cached AlphaFold entry for {uniprot_id}")
        with open(cache_path) as f:
            return json.load(f)

    url = ALPHAFOLD_API_URL.format(uniprot_id=uniprot_id)
    logger.info(f"Fetching AlphaFold entry for {uniprot_id}...")
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        logger.warning(f"AlphaFold fetch failed for {uniprot_id}: {e}")
        return None

    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Cached AlphaFold entry for {uniprot_id}")

    return data


def parse_alphafold_response(data) -> dict | None:
    """Return a dict with entry_id, pdb_url, uniprot from AlphaFold response.

    Returns None for empty/None input.
    """
    if not data:
        return None
    entry = data[0]
    return {
        "entry_id": entry.get("entryId"),
        "pdb_url": entry.get("pdbUrl"),
        "uniprot": entry.get("uniprotAccession"),
    }


def download_alphafold_pdb(entry: dict, gene_symbol: str) -> str | None:
    """Download an AlphaFold PDB file to local storage.

    Returns the local file path, or None on failure.
    """
    os.makedirs(STRUCTURES_DIR, exist_ok=True)
    entry_id = entry["entry_id"]
    dest = os.path.join(STRUCTURES_DIR, f"{gene_symbol}_{entry_id}.pdb")

    if os.path.exists(dest):
        logger.info(f"AlphaFold file already exists: {dest}")
        return dest

    url = entry["pdb_url"]
    logger.info(f"Downloading AlphaFold structure {entry_id} for {gene_symbol}...")
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException as e:
        logger.warning(f"Failed to download AlphaFold PDB {entry_id}: {e}")
        return None

    with open(dest, "w") as f:
        f.write(resp.text)
    logger.info(f"Saved AlphaFold PDB to {dest}")
    return dest


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def fetch_structure_for_gene(gene_symbol: str) -> dict | None:
    """Fetch the best available structure for a dockable gene.

    Priority: preferred_pdb → PDB search → AlphaFold fallback.

    Returns a dict with keys:
        gene_symbol, pdb_id, alphafold_id, structure_source,
        resolution_angstroms, file_path
    Returns None if gene is not in DOCKABLE_TARGETS.
    """
    entry = DOCKABLE_TARGETS.get(gene_symbol)
    if entry is None:
        logger.info(f"{gene_symbol} is not a dockable target — skipping")
        return None

    uniprot_id = entry["uniprot"]
    preferred_pdb = entry.get("preferred_pdb")

    # 1. Preferred PDB
    if preferred_pdb:
        logger.info(f"Using preferred PDB {preferred_pdb} for {gene_symbol}")
        file_path = download_pdb_file(preferred_pdb, gene_symbol)
        if file_path:
            return {
                "gene_symbol": gene_symbol,
                "pdb_id": preferred_pdb,
                "alphafold_id": None,
                "structure_source": "pdb_preferred",
                "resolution_angstroms": None,  # could be enriched later
                "file_path": file_path,
            }

    # 2. PDB search
    pdb_ids = search_pdb_for_gene(gene_symbol, uniprot_id)
    if pdb_ids:
        best_pdb = pdb_ids[0]
        logger.info(f"Using PDB search result {best_pdb} for {gene_symbol}")
        file_path = download_pdb_file(best_pdb, gene_symbol)
        if file_path:
            return {
                "gene_symbol": gene_symbol,
                "pdb_id": best_pdb,
                "alphafold_id": None,
                "structure_source": "pdb_search",
                "resolution_angstroms": None,
                "file_path": file_path,
            }

    # 3. AlphaFold fallback
    logger.info(f"Falling back to AlphaFold for {gene_symbol}")
    af_data = fetch_alphafold_entry(uniprot_id)
    af_entry = parse_alphafold_response(af_data)
    if af_entry:
        file_path = download_alphafold_pdb(af_entry, gene_symbol)
        if file_path:
            return {
                "gene_symbol": gene_symbol,
                "pdb_id": None,
                "alphafold_id": af_entry["entry_id"],
                "structure_source": "alphafold",
                "resolution_angstroms": None,
                "file_path": file_path,
            }

    logger.warning(f"No structure found for {gene_symbol}")
    return None


# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------

def store_protein_structure(conn, structure_data: dict) -> None:
    """Upsert a structure record into protein_structures table.

    Deletes any existing row for the gene_symbol before inserting, so only
    one row per gene is kept (the protein_structures table has no UNIQUE
    constraint on gene_symbol at the DDL level).
    """
    gene = structure_data["gene_symbol"]
    conn.execute("DELETE FROM protein_structures WHERE gene_symbol = ?", (gene,))
    conn.execute(
        """
        INSERT INTO protein_structures
            (gene_symbol, pdb_id, alphafold_id, structure_source,
             resolution_angstroms, file_path)
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        (
            gene,
            structure_data.get("pdb_id"),
            structure_data.get("alphafold_id"),
            structure_data["structure_source"],
            structure_data.get("resolution_angstroms"),
            structure_data["file_path"],
        ),
    )
    conn.commit()
    logger.info(f"Stored structure for {gene}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    conn = get_connection()
    init_db(conn)

    for gene_symbol in DOCKABLE_TARGETS:
        structure = fetch_structure_for_gene(gene_symbol)
        if structure:
            store_protein_structure(conn, structure)
        else:
            logger.warning(f"No structure stored for {gene_symbol}")

    conn.close()
    logger.info("Protein structure fetch complete.")


if __name__ == "__main__":
    main()
