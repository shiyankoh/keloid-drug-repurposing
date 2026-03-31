# Keloid Drug Repurposing Pipeline

A personal research pipeline to find FDA-approved drugs that might work for keloids by systematically analyzing molecular pathway overlaps.

**This is a hypothesis generator, not medical advice.** Discuss any candidate with a dermatologist before use.

## What it does

Keloids are driven by well-characterized molecular pathways (TGF-β, mTOR, PI3K/AKT, IL-6, collagen synthesis) that overlap with pathways targeted by drugs approved for other diseases. This pipeline surfaces those overlaps and ranks candidates by evidence strength.

1. **Seeds keloid biology** — 22 manually curated molecular targets across 10 pathways, enriched with data from [OpenTargets](https://www.opentargets.org/)
2. **Finds drug overlaps** — queries [DGIdb](https://www.dgidb.org/) for FDA-approved drugs that interact with keloid-associated genes
3. **Scores and ranks** — composite score based on pathway overlap count, evidence strength, and multi-target bonus
4. **Classifies safety** — green/yellow/red tiers based on whether a drug's side effect profile is appropriate for a benign condition
5. **Cross-references clinical trials** — checks [ClinicalTrials.gov](https://clinicaltrials.gov/) for existing keloid/hypertrophic scar trials

Output: a ranked Markdown report with 1,230 scored candidates, 12 actionable drugs, and 6 with existing clinical trials.

## Quick start

```bash
pip install -r requirements.txt

# Run the full pipeline
./run_pipeline.sh

# Or run individual steps
python3 -m scripts.01_seed_keloid_targets
python3 -m scripts.02_ingest_drugs
python3 -m scripts.03_score_candidates
python3 -m scripts.05_clinical_trials
python3 -m scripts.04_generate_report

# Run tests
python3 -m pytest tests/ -v
```

No API keys needed. All data sources are public and unauthenticated. API responses are cached locally so you can re-run without hitting external services.

## Tech stack

- **Python 3** — data processing
- **SQLite** — local database, zero config
- **pytest** — 57 unit tests
- **APIs:** OpenTargets (GraphQL), DGIdb (GraphQL), ClinicalTrials.gov (REST)

## Scoring formula

```
v1_score = (pathway_overlap × 0.40) + (evidence_strength × 0.40) + (multi_target_bonus × 0.20)
```

- **Pathway overlap** — how many curated keloid pathways does this drug target? (normalized 0–1)
- **Evidence strength** — average evidence score across matched gene targets (0–1)
- **Multi-target bonus** — 1.0 if the drug hits 3+ distinct curated pathways

## Severity tiers

Not every drug that hits keloid pathways is appropriate for treating a benign condition.

| Tier | Meaning | Example |
|------|---------|---------|
| 🟢 Green | Low risk, feasible for keloid use | Celecoxib, Dexamethasone, Aspirin |
| 🟡 Yellow | Moderate risk, needs careful evaluation | Imatinib, Rituximab, Lenalidomide |
| 🔴 Red | Oncology-grade toxicity, not appropriate | Sorafenib, Erlotinib |

## Known limitations

- Can't distinguish inhibition from activation — a drug that "touches" a pathway may not therapeutically modulate it in the right direction
- Scoring rewards breadth over precision — a highly targeted drug (e.g., verapamil) can be effective through a single mechanism but will score low
- DGIdb interaction types are often "unknown" — the nature of the drug-gene interaction isn't always specified
- Clinical trial *existence* is captured, but trial *results* are not analyzed
- Drug delivery (topical vs. injection vs. systemic) isn't modeled, which matters a lot for keloids

## Background

This pipeline is loosely modeled on [David Fajgenbaum's computational pharmacophenomics approach](https://everycure.org) for identifying repurposed drugs for rare diseases.

## Feedback welcome

If you're a dermatologist, computational biologist, or anyone who works on keloids — I'd love your thoughts on the approach, the candidate list, or the experiment design. Reach out at keloid-research@agentmail.to.

## License

MIT
