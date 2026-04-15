"""Import directionality annotations from ralph's JSON output into drug_targets."""

import json
import logging
import os

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

ANNOTATIONS_PATH = os.path.join(
    os.path.dirname(__file__), "..", "data", "directionality_annotations.json"
)


def import_annotations(conn, annotations_path=None):
    """Read annotation JSON and upsert action_direction into drug_targets.

    Returns the number of rows updated.
    """
    path = annotations_path or ANNOTATIONS_PATH
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Annotations file not found: {path}. Run the ralph loop first."
        )

    with open(path) as f:
        annotations = json.load(f)

    REQUIRED_FIELDS = {"drug_name", "gene_symbol", "action_direction"}
    for i, ann in enumerate(annotations):
        missing = REQUIRED_FIELDS - ann.keys()
        if missing:
            raise ValueError(
                f"Annotation entry {i} missing fields: {missing}. "
                f"Check if the ralph loop completed successfully."
            )

    updated = 0
    for ann in annotations:
        drug_name = ann["drug_name"]
        gene_symbol = ann["gene_symbol"]
        action_direction = ann["action_direction"]

        cursor = conn.execute(
            """
            UPDATE drug_targets
            SET action_direction = ?
            WHERE gene_symbol = ?
              AND drug_id IN (SELECT id FROM drugs WHERE drug_name = ?)
            """,
            (action_direction, gene_symbol, drug_name),
        )

        if cursor.rowcount > 0:
            updated += cursor.rowcount
            logger.info(f"  {drug_name} → {gene_symbol}: {action_direction}")
        else:
            logger.warning(f"  No match for {drug_name} / {gene_symbol} — skipping")

    conn.commit()
    logger.info(f"Imported {updated} directionality annotations")
    return updated


def main():
    from scripts.db import get_connection, migrate_directionality_columns

    conn = get_connection()
    migrate_directionality_columns(conn)
    import_annotations(conn)
    conn.close()
