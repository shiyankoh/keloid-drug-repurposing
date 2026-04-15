# SIROLIMUS — Directionality Research

## Drug Class
mTOR inhibitor (macrolide immunosuppressant / rapamycin)

## ChEMBL ID
CHEMBL413

## Mechanism Overview
Sirolimus (rapamycin) binds FKBP12 (FK506-binding protein 1A); the FKBP12-sirolimus complex directly inhibits mTORC1 (mechanistic target of rapamycin complex 1). ChEMBL lists the primary mechanism as "FK506-binding protein 1A inhibitor" (CHEMBL1902). mTOR is the functional pharmacological target — sirolimus is the canonical mTOR inhibitor and the founding member of the "rapalog" class.

## Gene-by-Gene Directionality

| Gene | Direction | Confidence | Net Effect | Rationale |
|------|-----------|------------|------------|-----------|
| MTOR | inhibitor | high | beneficial | Canonical target. FKBP12-sirolimus complex directly inhibits mTORC1 kinase activity, blocking S6K1 and 4E-BP1 phosphorylation. Well-validated in fibrosis models — reduces fibroblast proliferation and collagen synthesis. |
| PIK3CA | modulator | medium | unclear | Complex relationship: mTOR is downstream of PI3K, but mTORC1 inhibition relieves S6K-mediated negative feedback on IRS-1, potentially activating PI3K/AKT (O'Reilly et al. 2006). Net PI3K pathway activity depends on cell context. In fibrosis models, rapamycin generally suppresses proliferation despite feedback. |
| EGFR | inhibitor | low | beneficial | Indirect effect — mTOR inhibition reduces cap-dependent translation, which can lower EGFR protein levels. Some evidence that rapamycin sensitizes cells to EGFR inhibitors. No direct EGFR kinase inhibition. |
| JAK2 | inhibitor | low | beneficial | Indirect — mTOR inhibition has been shown to reduce JAK-STAT signaling output in myeloproliferative models. Rapamycin suppresses JAK2-V617F-driven proliferation in synergy with JAK2 inhibitors. No direct JAK2 kinase inhibition. |

## Summary
Sirolimus is the canonical mTOR inhibitor with high-confidence direct inhibition of MTOR (mTORC1). Its effects on the other three keloid targets are indirect: PIK3CA is complicated by a well-documented feedback activation loop (mTORC1 inhibition → IRS-1 derepression → PI3K activation), making the net direction unclear. EGFR and JAK2 effects are downstream consequences of reduced mTOR signaling with low confidence. Overall, sirolimus has a favorable profile for keloid treatment given its strong anti-proliferative and anti-fibrotic properties.

## Sources
- ChEMBL: CHEMBL413 (FK506-binding protein 1A inhibitor)
- PubMed:12461523 — Rapamycin inhibition of mTOR (Sarbassov et al., molecular biology review)
- PubMed:16489089 — O'Reilly et al., mTOR inhibition induces upstream PI3K/AKT activation via feedback loop
- PubMed:20068183 — Rapamycin in fibrosis models: reduces fibroblast proliferation and collagen
- PubMed:22031083 — mTOR pathway in fibroproliferative disorders
- PubMed:20829794 — Rapamycin synergy with JAK2 inhibitors in myeloproliferative neoplasms
- PubMed:18280086 — mTOR inhibition and EGFR signaling crosstalk
