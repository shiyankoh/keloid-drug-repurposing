"""Score drug candidates by keloid pathway overlap."""

import logging
from scripts.db import get_connection, init_db

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

W_OVERLAP = 0.40
W_EVIDENCE = 0.40
W_MULTITARGET = 0.20
MULTI_TARGET_THRESHOLD = 3
UNCURATED_PATHWAY = "OpenTargets_association"


def query_overlaps(conn):
    """Find all drugs that target keloid-associated genes.
    Returns one row per drug-gene match with pathway and evidence info."""
    cursor = conn.execute("""
        SELECT d.id as drug_id, d.drug_name, d.generic_name,
               d.original_indication, d.mechanism_of_action,
               dt.gene_symbol, dt.action_type,
               kt.pathway, kt.evidence_strength
        FROM drugs d
        JOIN drug_targets dt ON d.id = dt.drug_id
        JOIN keloid_targets kt ON dt.gene_symbol = kt.gene_symbol
        WHERE d.approval_status = 'approved'
        ORDER BY d.drug_name, kt.pathway
    """)
    return [dict(row) for row in cursor.fetchall()]


def normalize_overlap_counts(scores):
    """Normalize pathway_overlap_count to 0.0-1.0 range."""
    if not scores:
        return []
    max_count = max(s["pathway_overlap_count"] for s in scores)
    for s in scores:
        s["normalized_overlap"] = s["pathway_overlap_count"] / max_count
    return scores


def compute_scores(overlaps):
    """Aggregate overlaps into per-drug scores."""
    drug_data = {}
    for row in overlaps:
        name = row["drug_name"]
        if name not in drug_data:
            drug_data[name] = {
                "drug_name": name,
                "generic_name": row["generic_name"],
                "original_indication": row["original_indication"],
                "mechanism_of_action": row["mechanism_of_action"],
                "pathways": set(),
                "genes": [],
                "evidence_values": [],
                "gene_details": [],
            }
        drug_data[name]["pathways"].add(row["pathway"])
        drug_data[name]["genes"].append(row["gene_symbol"])
        drug_data[name]["evidence_values"].append(row["evidence_strength"])
        drug_data[name]["gene_details"].append({
            "gene_symbol": row["gene_symbol"],
            "pathway": row["pathway"],
            "action_type": row["action_type"],
            "evidence_strength": row["evidence_strength"],
        })

    scores = []
    for name, data in drug_data.items():
        curated_pathways = {p for p in data["pathways"] if p != UNCURATED_PATHWAY}
        pathway_count = len(curated_pathways)
        avg_evidence = sum(data["evidence_values"]) / len(data["evidence_values"])
        multi_bonus = 1.0 if pathway_count >= MULTI_TARGET_THRESHOLD else 0.0

        scores.append({
            "drug_name": name,
            "generic_name": data["generic_name"],
            "original_indication": data["original_indication"],
            "mechanism_of_action": data["mechanism_of_action"],
            "pathway_overlap_count": pathway_count,
            "pathways": sorted(data["pathways"]),
            "avg_evidence_strength": round(avg_evidence, 4),
            "multi_target_bonus": multi_bonus,
            "gene_details": data["gene_details"],
            "composite_score": 0.0,
        })

    scores = normalize_overlap_counts(scores)
    for s in scores:
        s["composite_score"] = round(
            s["normalized_overlap"] * W_OVERLAP
            + s["avg_evidence_strength"] * W_EVIDENCE
            + s["multi_target_bonus"] * W_MULTITARGET,
            4,
        )

    scores.sort(key=lambda s: s["composite_score"], reverse=True)
    return scores


def main():
    conn = get_connection()

    overlaps = query_overlaps(conn)
    if not overlaps:
        logger.warning("No drug-keloid overlaps found. Check that steps 01 and 02 ran successfully.")
        conn.close()
        return

    scores = compute_scores(overlaps)
    logger.info(f"Scored {len(scores)} drug candidates")

    print(f"\n{'Rank':<6}{'Drug':<25}{'Score':<8}{'Pathways':<10}{'Avg Evidence':<14}{'Multi-target'}")
    print("-" * 75)
    for i, s in enumerate(scores[:10], 1):
        print(
            f"{i:<6}{s['drug_name']:<25}{s['composite_score']:<8.4f}"
            f"{s['pathway_overlap_count']:<10}{s['avg_evidence_strength']:<14.4f}"
            f"{'YES' if s['multi_target_bonus'] else 'no'}"
        )

    import json
    import os
    output_dir = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
    os.makedirs(output_dir, exist_ok=True)
    scores_path = os.path.join(output_dir, "scores.json")
    with open(scores_path, "w") as f:
        json.dump(scores, f, indent=2, default=list)
    logger.info(f"Full scores saved to {scores_path}")

    conn.close()


if __name__ == "__main__":
    main()
