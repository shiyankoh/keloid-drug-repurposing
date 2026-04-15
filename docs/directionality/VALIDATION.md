# Directionality Research Validation

## Summary

| Metric | Expected | Actual | Status |
|--------|----------|--------|--------|
| Research docs in `docs/directionality/` | 12 | 12 | PASS |
| Entries in `directionality_annotations.json` | 38 | 38 | PASS |
| All entries have 8 required fields | Yes | Yes | PASS |
| Valid `action_direction` values | Yes | Yes | PASS |
| Valid `confidence` values | Yes | Yes | PASS |
| Valid `net_effect` values | Yes | Yes | PASS |
| TNF entries have `net_effect=unclear` | Yes | Yes | PASS |

## Drug-Gene Pair Counts

| Drug | Pairs | Expected |
|------|-------|----------|
| LENALIDOMIDE | 4 | 4 (TNF, VEGFA, KDR, CTNNB1) |
| RITUXIMAB | 3 | 3 (IL6, JAK2, TGFB1) |
| IMATINIB | 4 | 4 (EGFR, JAK2, SMAD3, CTNNB1) |
| CELECOXIB | 4 | 4 (STAT3, PIK3CA, VEGFA, CTNNB1) |
| PEGINTERFERON ALFA-2B | 3 | 3 (EGFR, IL6, JAK2) |
| DEXAMETHASONE | 4 | 4 (SMAD3, VEGFA, KDR, CTNNB1) |
| THALIDOMIDE | 3 | 3 (TNF, VEGFA, CTNNB1) |
| SIROLIMUS | 4 | 4 (MTOR, PIK3CA, EGFR, JAK2) |
| ASPIRIN | 3 | 3 (PIK3CA, TGFB1, CYP24A1) |
| CUCURBITACIN B | 2 | 2 (STAT3, SMAD3) |
| MILTEFOSINE | 2 | 2 (STAT3, SMAD3) |
| TRIFLUOPERAZINE | 2 | 2 (STAT3, SMAD3) |
| **Total** | **38** | **38** |

## Action Direction Distribution

| Direction | Count |
|-----------|-------|
| inhibitor | 31 |
| activator | 2 |
| modulator | 1 |
| unknown | 4 |

## Unknown Pairs Flagged

These 4 drug-gene pairs could not be assigned a clear direction:

1. **RITUXIMAB-JAK2** — No direct or well-documented indirect mechanism (low confidence)
2. **IMATINIB-EGFR** — Imatinib does not inhibit EGFR; distinct kinase family (low confidence)
3. **ASPIRIN-CYP24A1** — No published evidence of aspirin modulating CYP24A1 (low confidence)
4. **TRIFLUOPERAZINE-STAT3** — Conflicting theoretical vs experimental evidence (low confidence)

## Unfavorable Net Effects

2 entries where the drug activates a pro-keloid target:

1. **PEGINTERFERON ALFA-2B-JAK2** — Type I IFN activates JAK2 via IFNAR signaling
2. **MILTEFOSINE-STAT3** — Compensatory JAK2/STAT3 activation via Akt inhibition

## Confidence Distribution

| Confidence | Count |
|------------|-------|
| high | 5 |
| medium | 21 |
| low | 12 |

## Research Docs

All 12 files present in `docs/directionality/`:
ASPIRIN.md, CELECOXIB.md, CUCURBITACIN_B.md, DEXAMETHASONE.md, IMATINIB.md, LENALIDOMIDE.md, MILTEFOSINE.md, PEGINTERFERON_ALFA-2B.md, RITUXIMAB.md, SIROLIMUS.md, THALIDOMIDE.md, TRIFLUOPERAZINE.md

## Completeness

All 12 drugs researched. All 38 drug-gene pairs annotated. All required fields populated with valid values. Dataset is ready for import by `scripts/08_import_directionality.py`.
