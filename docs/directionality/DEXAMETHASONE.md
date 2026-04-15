# DEXAMETHASONE — Directionality Research

## Drug Class
Corticosteroid / Glucocorticoid receptor agonist

## ChEMBL ID
CHEMBL384467

## Mechanism of Action
Dexamethasone is a synthetic glucocorticoid that activates the glucocorticoid receptor (GR/NR3C1). GR translocates to the nucleus and modulates transcription — suppressing pro-inflammatory genes (via NF-κB and AP-1 inhibition) and inducing anti-inflammatory genes. It is the most potent commonly used corticosteroid.

ChEMBL lists one mechanism: **Glucocorticoid receptor agonist** (target CHEMBL2034). GR is not in our 22-gene keloid target list, so all effects on SMAD3, VEGFA, KDR, and CTNNB1 are indirect, mediated through GR-dependent transcriptional regulation.

## Gene-by-Gene Directionality

| Gene | Direction | Confidence | Net Effect | Source |
|------|-----------|------------|------------|--------|
| SMAD3 | inhibitor | medium | beneficial | PubMed:15681343, PubMed:11483860 |
| VEGFA | inhibitor | high | beneficial | PubMed:12042317, PubMed:16081860 |
| KDR | inhibitor | medium | beneficial | PubMed:15520188, PubMed:17035535 |
| CTNNB1 | inhibitor | low | beneficial | PubMed:15878872, PubMed:20133788 |

### SMAD3
Dexamethasone inhibits TGF-β/SMAD signaling in fibroblasts. GR activation interferes with SMAD3-dependent transcription through direct GR-SMAD3 protein interaction and by inducing SMAD7 (an inhibitory SMAD). Demonstrated in lung fibroblasts and renal mesangial cells — dexamethasone reduces TGF-β1-induced SMAD3 phosphorylation and nuclear translocation. Medium confidence because the mechanism is indirect but well-characterized in fibrotic contexts.

### VEGFA
Dexamethasone is a potent inhibitor of VEGF expression. This is clinically exploited — intravitreal dexamethasone implants (Ozurdex) treat macular edema partly through VEGF suppression. GR activation suppresses VEGF transcription via NF-κB and AP-1 inhibition, and reduces HIF-1α-mediated VEGF induction under hypoxia. High confidence due to extensive clinical and preclinical evidence.

### KDR
Corticosteroids including dexamethasone downregulate VEGFR2/KDR expression on endothelial cells. This contributes to their anti-angiogenic effects beyond VEGF suppression alone. Demonstrated in dermal microvascular endothelial cells and retinal models. Medium confidence — well-documented but most studies are in vascular rather than fibroblast contexts.

### CTNNB1
Dexamethasone modulates Wnt/β-catenin signaling, generally suppressing it. GR activation can induce GSK-3β activity (promoting β-catenin degradation) and directly interact with TCF/LEF transcription factors to inhibit β-catenin-dependent gene expression. However, the literature is mixed — in osteoblasts, GR suppresses Wnt signaling (contributing to glucocorticoid-induced osteoporosis), while in some cancer contexts effects vary. Low confidence for keloid fibroblasts specifically, though the anti-proliferative direction is consistent.

## Summary
Dexamethasone acts through the glucocorticoid receptor to broadly suppress pro-inflammatory and pro-fibrotic signaling. All four keloid-relevant targets (SMAD3, VEGFA, KDR, CTNNB1) are inhibited, yielding a uniformly beneficial profile. VEGFA inhibition has the strongest evidence (clinical use in ophthalmology), while CTNNB1 inhibition is least certain due to context-dependent Wnt/GR crosstalk.

## Sources
- ChEMBL: CHEMBL384467 (Glucocorticoid receptor agonist)
- PubMed:15681343 — GR-SMAD3 interaction in fibroblasts
- PubMed:11483860 — Dexamethasone induces SMAD7, inhibits TGF-β signaling
- PubMed:12042317 — Dexamethasone suppresses VEGF in retinal pigment epithelium
- PubMed:16081860 — Corticosteroid regulation of VEGF expression
- PubMed:15520188 — Glucocorticoid effects on VEGFR2/KDR expression
- PubMed:17035535 — Anti-angiogenic mechanisms of corticosteroids
- PubMed:15878872 — GR-mediated suppression of Wnt/β-catenin signaling
- PubMed:20133788 — Glucocorticoids and GSK-3β/β-catenin in bone
