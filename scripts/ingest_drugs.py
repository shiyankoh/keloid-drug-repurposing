"""Ingest drug-gene interactions from DGIdb for all keloid target genes."""

import json
import logging
import os
import requests
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DGIDB_URL = "https://dgidb.org/api/graphql"
RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw")

DGIDB_QUERY = """
query drugInteractions($genes: [String!]!) {
  genes(names: $genes) {
    nodes {
      name
      longName
      interactions {
        drug {
          name
          approved
          conceptId
        }
        interactionScore
        interactionTypes {
          type
          directionality
        }
      }
    }
  }
}
"""


def fetch_dgidb(gene_symbols, cache_path=None):
    """Fetch drug-gene interactions from DGIdb. Uses cache if available."""
    cache = cache_path or os.path.join(RAW_DIR, "dgidb_interactions.json")

    if os.path.exists(cache):
        logger.info(f"Loading cached DGIdb response from {cache}")
        with open(cache) as f:
            return json.load(f)

    logger.info(f"Fetching interactions for {len(gene_symbols)} genes from DGIdb...")
    resp = requests.post(
        DGIDB_URL,
        json={"query": DGIDB_QUERY, "variables": {"genes": gene_symbols}},
        timeout=60,
    )
    resp.raise_for_status()
    data = resp.json()

    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Cached DGIdb response to {cache}")

    return data


def parse_dgidb_response(response_json):
    """Parse DGIdb response into deduplicated drugs list and interactions list.
    Filters to FDA-approved drugs only."""
    nodes = response_json.get("data", {}).get("genes", {}).get("nodes", [])

    seen_drugs = {}
    interactions = []

    for gene_node in nodes:
        gene_symbol = gene_node["name"]
        for interaction in gene_node.get("interactions", []):
            drug = interaction["drug"]
            if not drug.get("approved", False):
                continue

            drug_name = drug["name"]
            drug_key = drug_name.lower()

            if drug_key not in seen_drugs:
                seen_drugs[drug_key] = {
                    "drug_name": drug_name,
                    "generic_name": drug_name.lower(),
                    "approval_status": "approved",
                    "original_indication": None,
                    "mechanism_of_action": None,
                }

            action_types = interaction.get("interactionTypes", [])
            action_type = action_types[0]["type"] if action_types else "unknown"

            interactions.append({
                "drug_name": drug_name,
                "gene_symbol": gene_symbol,
                "action_type": action_type,
                "source": "DGIdb",
            })

    drugs = list(seen_drugs.values())
    return drugs, interactions


def insert_drugs_and_targets(conn, drugs, interactions):
    """Insert drugs and drug_targets into SQLite."""
    drug_id_map = {}
    for d in drugs:
        cursor = conn.execute(
            """INSERT OR IGNORE INTO drugs (drug_name, generic_name, approval_status,
               original_indication, mechanism_of_action)
               VALUES (?, ?, ?, ?, ?)""",
            (d["drug_name"], d["generic_name"], d["approval_status"],
             d["original_indication"], d["mechanism_of_action"]),
        )
        if cursor.lastrowid:
            drug_id_map[d["drug_name"]] = cursor.lastrowid
        else:
            row = conn.execute(
                "SELECT id FROM drugs WHERE generic_name = ?", (d["generic_name"],)
            ).fetchone()
            drug_id_map[d["drug_name"]] = row["id"]

    for i in interactions:
        drug_id = drug_id_map[i["drug_name"]]
        conn.execute(
            """INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source)
               VALUES (?, ?, ?, ?)""",
            (drug_id, i["gene_symbol"], i["action_type"], i["source"]),
        )

    conn.commit()
    drug_count = conn.execute("SELECT COUNT(*) as cnt FROM drugs").fetchone()["cnt"]
    target_count = conn.execute("SELECT COUNT(*) as cnt FROM drug_targets").fetchone()["cnt"]
    logger.info(f"Inserted {drug_count} drugs with {target_count} drug-target interactions")


def main():
    conn = get_connection()
    init_db(conn)

    conn.execute("DELETE FROM drug_targets")
    conn.execute("DELETE FROM drugs")
    conn.commit()

    cursor = conn.execute("SELECT gene_symbol FROM keloid_targets")
    gene_symbols = [row["gene_symbol"] for row in cursor.fetchall()]

    if not gene_symbols:
        logger.warning("No keloid targets found. Run 01_seed_keloid_targets.py first.")
        conn.close()
        return

    logger.info(f"Querying DGIdb for {len(gene_symbols)} keloid target genes")

    response = fetch_dgidb(gene_symbols)
    drugs, interactions = parse_dgidb_response(response)
    logger.info(f"Found {len(drugs)} approved drugs with {len(interactions)} interactions")

    insert_drugs_and_targets(conn, drugs, interactions)
    conn.close()


if __name__ == "__main__":
    main()
