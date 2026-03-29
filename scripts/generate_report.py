"""Generate a Markdown report from scored drug candidates."""

import json
import logging
import os
from datetime import date

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

SCORES_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "scores.json")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output")

TIER_LABELS = {
    "green": "Low risk",
    "yellow": "Moderate risk",
    "red": "High risk",
    "unclassified": "Not yet classified",
}

TIER_EMOJI = {
    "green": "G",
    "yellow": "Y",
    "red": "R",
    "unclassified": "?",
}


def generate_report(scores, top_n=20):
    """Generate Markdown report string from scored candidates."""
    today = date.today().isoformat()

    # Separate actionable candidates (green/yellow with 2+ curated pathways)
    actionable = [
        s for s in scores
        if s.get("severity_tier") in ("green", "yellow")
        and s["pathway_overlap_count"] >= 2
    ]

    lines = [
        f"# Keloid Drug Repurposing Candidates",
        f"",
        f"*Generated {today} — {len(scores)} candidates scored*",
        f"",
        f"**Scoring formula:** `v1_score = overlap x 0.40 + evidence x 0.40 + multi_target x 0.20`",
        f"",
        f"**Severity tiers:** [G] Green = low risk, feasible for keloid use | "
        f"[Y] Yellow = moderate risk, needs careful evaluation | "
        f"[R] Red = high risk, oncology-grade toxicity | "
        f"[?] Unclassified",
        f"",
        f"> This is a hypothesis generator, not medical advice. Discuss any candidate with a dermatologist before use.",
        f"",
        f"---",
        f"",
    ]

    # Actionable candidates section (the most useful part)
    lines.extend([
        f"## Actionable Candidates ({len(actionable)} drugs)",
        f"",
        f"These drugs have green/yellow safety profiles AND hit 2+ curated keloid pathways.",
        f"Sorted by composite score.",
        f"",
        f"| Rank | Drug | Score | Pathways | Avg Evidence | Tier | Safety Notes |",
        f"|------|------|-------|----------|-------------|------|-------------|",
    ])

    for i, s in enumerate(actionable, 1):
        tier = s.get("severity_tier", "unclassified")
        tier_tag = f"[{TIER_EMOJI[tier]}]"
        notes = s.get("severity_notes", "")
        # Truncate notes for table
        short_notes = notes[:60] + "..." if len(notes) > 60 else notes
        lines.append(
            f"| {i} | {s['drug_name']} | {s['composite_score']:.4f} | "
            f"{s['pathway_overlap_count']} | {s['avg_evidence_strength']:.4f} | "
            f"{tier_tag} | {short_notes} |"
        )

    lines.extend(["", "---", ""])

    # Clinical trials section
    drugs_with_trials = [s for s in scores if s.get("clinical_trials")]
    if drugs_with_trials:
        lines.extend([
            f"## Clinical Trial Evidence ({len(drugs_with_trials)} drugs with keloid/scar trials)",
            f"",
            f"These candidates have existing or completed trials on ClinicalTrials.gov for keloid or hypertrophic scar.",
            f"",
            f"| Drug | Tier | Score | # Trials | Phases | Status |",
            f"|------|------|-------|----------|--------|--------|",
        ])
        for s in drugs_with_trials:
            trials = s["clinical_trials"]
            tier_tag = f"[{TIER_EMOJI.get(s.get('severity_tier', 'unclassified'), '?')}]"
            phases = ", ".join(sorted(set(t["phase"] for t in trials)))
            statuses = ", ".join(sorted(set(t["status"] for t in trials)))
            lines.append(
                f"| {s['drug_name']} | {tier_tag} | {s['composite_score']:.4f} | "
                f"{len(trials)} | {phases} | {statuses} |"
            )
        lines.extend([""])

        # List individual trials
        for s in drugs_with_trials:
            lines.append(f"**{s['drug_name']}:**")
            for t in s["clinical_trials"]:
                status_str = t['status'].replace('_', ' ').title()
                lines.append(
                    f"- [{t['nct_id']}](https://clinicaltrials.gov/study/{t['nct_id']}) — "
                    f"{t['title']} ({t['phase']}, {status_str})"
                )
            lines.append("")

        lines.extend(["---", ""])

    # Full top N (all tiers)
    lines.extend([
        f"## All Top {min(top_n, len(scores))} Candidates (any tier)",
        f"",
        f"| Rank | Drug | Score | Pathways | Avg Evidence | Multi-target | Tier |",
        f"|------|------|-------|----------|-------------|-------------|------|",
    ])

    for i, s in enumerate(scores[:top_n], 1):
        mt = "Yes" if s["multi_target_bonus"] else "-"
        tier = s.get("severity_tier", "unclassified")
        tier_tag = f"[{TIER_EMOJI[tier]}]"
        lines.append(
            f"| {i} | {s['drug_name']} | {s['composite_score']:.4f} | "
            f"{s['pathway_overlap_count']} | {s['avg_evidence_strength']:.4f} | "
            f"{mt} | {tier_tag} |"
        )

    lines.extend(["", "---", ""])

    # Detailed briefs for top actionable candidates
    brief_candidates = actionable[:10] if actionable else scores[:10]
    brief_label = "Top Actionable" if actionable else "Top 10"
    lines.extend([f"## Detailed Briefs ({brief_label})", ""])

    for i, s in enumerate(brief_candidates, 1):
        tier = s.get("severity_tier", "unclassified")
        lines.append(f"### {i}. {s['drug_name']} [{TIER_EMOJI[tier]}]")
        lines.append("")
        if s.get("severity_notes"):
            lines.append(f"**Safety:** {s['severity_notes']}")
        if s.get("original_indication"):
            lines.append(f"**Approved for:** {s['original_indication']}")
        if s.get("mechanism_of_action"):
            lines.append(f"**Mechanism:** {s['mechanism_of_action']}")
        lines.append(f"**Composite score:** {s['composite_score']:.4f}")
        lines.append(f"**Distinct keloid pathways:** {s['pathway_overlap_count']} ({', '.join(s['pathways'])})")
        lines.append("")

        # Only show curated pathway gene details (skip OpenTargets noise in briefs)
        curated_details = [g for g in s.get("gene_details", []) if g["pathway"] != "OpenTargets_association"]
        ot_details = [g for g in s.get("gene_details", []) if g["pathway"] == "OpenTargets_association"]

        if curated_details:
            lines.append("| Gene | Pathway | Action | Evidence |")
            lines.append("|------|---------|--------|----------|")
            for g in curated_details:
                lines.append(
                    f"| {g['gene_symbol']} | {g['pathway']} | "
                    f"{g['action_type']} | {g['evidence_strength']:.2f} |"
                )
            if ot_details:
                lines.append(f"")
                lines.append(f"*Plus {len(ot_details)} additional OpenTargets gene associations (low confidence, omitted for clarity)*")

        # Show trial info if available
        drug_trials = s.get("clinical_trials", [])
        if drug_trials:
            lines.append("")
            lines.append(f"**Clinical trials ({len(drug_trials)}):**")
            for t in drug_trials:
                status_str = t['status'].replace('_', ' ').title()
                lines.append(
                    f"- [{t['nct_id']}](https://clinicaltrials.gov/study/{t['nct_id']}) — "
                    f"{t['title']} ({t['phase']}, {status_str})"
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
