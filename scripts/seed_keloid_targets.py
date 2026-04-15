"""Seed keloid_targets table from curated CSV + OpenTargets API."""

import csv
import json
import logging
import os
import requests
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
KELOID_DISEASE_ID = "EFO_0004212"
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
SEED_CSV = os.path.join(DATA_DIR, "seed_targets.csv")
PAGE_SIZE = 500

OPENTARGETS_QUERY = """
query keloidTargets($diseaseId: String!, $size: Int!) {
  disease(efoId: $diseaseId) {
    id
    name
    associatedTargets(page: {index: 0, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""


def parse_seed_csv(file_obj):
    """Parse the seed CSV into a list of target dicts."""
    reader = csv.DictReader(file_obj)
    targets = []
    for row in reader:
        targets.append({
            "gene_symbol": row["gene_symbol"],
            "target_name": row["target_name"],
            "ensembl_id": None,
            "pathway": row["pathway"],
            "evidence_type": row["evidence_type"],
            "evidence_strength": float(row["evidence_strength"]),
            "source": row["source"],
            "target_role": row.get("target_role", "pro_keloid"),
        })
    return targets


def parse_opentargets_response(response_json):
    """Parse OpenTargets GraphQL response into a list of target dicts."""
    disease = response_json.get("data", {}).get("disease")
    if disease is None:
        logger.warning("OpenTargets returned no disease data for keloid query")
        return []

    assoc = disease.get("associatedTargets", {})
    rows = assoc.get("rows", [])
    if not rows:
        logger.warning("OpenTargets returned 0 associated targets for keloid")
        return []

    targets = []
    for row in rows:
        target = row["target"]
        targets.append({
            "gene_symbol": target["approvedSymbol"],
            "target_name": target["approvedName"],
            "ensembl_id": target["id"],
            "pathway": "OpenTargets_association",
            "evidence_type": "computed",
            "evidence_strength": round(row["score"], 4),
            "source": "opentargets",
        })
    return targets


def fetch_opentargets(cache_path=None):
    """Fetch keloid targets from OpenTargets. Uses cache if available."""
    cache = cache_path or os.path.join(RAW_DIR, "opentargets_keloid.json")

    if os.path.exists(cache):
        logger.info(f"Loading cached OpenTargets response from {cache}")
        with open(cache) as f:
            return json.load(f)

    logger.info("Fetching keloid targets from OpenTargets API...")
    resp = requests.post(
        OPENTARGETS_URL,
        json={
            "query": OPENTARGETS_QUERY,
            "variables": {"diseaseId": KELOID_DISEASE_ID, "size": PAGE_SIZE},
        },
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()

    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Cached OpenTargets response to {cache}")

    return data


def insert_targets(conn, csv_targets, api_targets):
    """Insert targets into keloid_targets. CSV seeds are inserted first
    so they win on gene_symbol UNIQUE conflict (INSERT OR IGNORE)."""
    csv_sql = """
        INSERT OR IGNORE INTO keloid_targets
        (gene_symbol, target_name, ensembl_id, pathway, evidence_type, evidence_strength, source, target_role)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """
    api_sql = """
        INSERT OR IGNORE INTO keloid_targets
        (gene_symbol, target_name, ensembl_id, pathway, evidence_type, evidence_strength, source)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    """
    for t in csv_targets:
        conn.execute(csv_sql, (
            t["gene_symbol"], t["target_name"], t["ensembl_id"],
            t["pathway"], t["evidence_type"], t["evidence_strength"], t["source"],
            t.get("target_role", "pro_keloid"),
        ))

    for t in api_targets:
        conn.execute(api_sql, (
            t["gene_symbol"], t["target_name"], t["ensembl_id"],
            t["pathway"], t["evidence_type"], t["evidence_strength"], t["source"],
        ))

    conn.commit()
    cursor = conn.execute("SELECT COUNT(*) as cnt FROM keloid_targets")
    count = cursor.fetchone()["cnt"]
    logger.info(f"Inserted {count} keloid targets")


def main():
    conn = get_connection()
    init_db(conn)

    conn.execute("DELETE FROM keloid_targets")
    conn.commit()

    with open(SEED_CSV) as f:
        csv_targets = parse_seed_csv(f)
    logger.info(f"Loaded {len(csv_targets)} targets from seed CSV")

    api_response = fetch_opentargets()
    api_targets = parse_opentargets_response(api_response)
    logger.info(f"Parsed {len(api_targets)} targets from OpenTargets")

    insert_targets(conn, csv_targets, api_targets)
    conn.close()


if __name__ == "__main__":
    main()
