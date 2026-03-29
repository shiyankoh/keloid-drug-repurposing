"""Generate a Markdown report from scored drug candidates."""

import json
import logging
import os
from datetime import date

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

SCORES_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "scores.json")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output")


def generate_report(scores, top_n=20):
    """Generate Markdown report string from scored candidates."""
    today = date.today().isoformat()
    lines = [
        f"# Keloid Drug Repurposing Candidates",
        f"",
        f"*Generated {today} — {len(scores)} candidates scored*",
        f"",
        f"**Scoring formula:** `v1_score = overlap×0.40 + evidence×0.40 + multi_target×0.20`",
        f"",
        f"> This is a hypothesis generator, not medical advice. Discuss any candidate with a dermatologist before use.",
        f"",
        f"---",
        f"",
        f"## Top {min(top_n, len(scores))} Candidates",
        f"",
        f"| Rank | Drug | Score | Pathways Hit | Avg Evidence | Multi-target |",
        f"|------|------|-------|-------------|-------------|-------------|",
    ]

    for i, s in enumerate(scores[:top_n], 1):
        mt = "Yes" if s["multi_target_bonus"] else "-"
        lines.append(
            f"| {i} | {s['drug_name']} | {s['composite_score']:.4f} | "
            f"{s['pathway_overlap_count']} | {s['avg_evidence_strength']:.4f} | {mt} |"
        )

    lines.extend(["", "---", ""])

    lines.extend([f"## Detailed Briefs (Top 10)", ""])
    for i, s in enumerate(scores[:10], 1):
        lines.append(f"### {i}. {s['drug_name']}")
        lines.append("")
        if s.get("original_indication"):
            lines.append(f"**Approved for:** {s['original_indication']}")
        if s.get("mechanism_of_action"):
            lines.append(f"**Mechanism:** {s['mechanism_of_action']}")
        lines.append(f"**Composite score:** {s['composite_score']:.4f}")
        lines.append(f"**Distinct keloid pathways:** {s['pathway_overlap_count']} ({', '.join(s['pathways'])})")
        lines.append("")

        lines.append("| Gene | Pathway | Action | Evidence |")
        lines.append("|------|---------|--------|----------|")
        for g in s.get("gene_details", []):
            lines.append(
                f"| {g['gene_symbol']} | {g['pathway']} | "
                f"{g['action_type']} | {g['evidence_strength']:.2f} |"
            )
        lines.extend([
            "",
            "**Next steps:** Discuss with dermatologist. Search PubMed for "
            f"`\"{s['drug_name'].lower()}\" AND (keloid OR fibrosis)`. "
            "Check ClinicalTrials.gov for existing trials.",
            "",
            "---",
            "",
        ])

    return "\n".join(lines)


def main():
    if not os.path.exists(SCORES_PATH):
        logger.error(f"Scores file not found at {SCORES_PATH}. Run 03_score_candidates.py first.")
        return

    with open(SCORES_PATH) as f:
        scores = json.load(f)

    if not scores:
        logger.warning("No scored candidates found.")
        return

    report = generate_report(scores)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_path = os.path.join(OUTPUT_DIR, "report.md")
    with open(output_path, "w") as f:
        f.write(report)

    logger.info(f"Report written to {output_path} ({len(scores)} candidates)")


if __name__ == "__main__":
    main()
