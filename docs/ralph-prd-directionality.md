# Keloid Directionality Research Agent

## Goal

For each of the 12 actionable keloid drug candidates below, determine whether each drug INHIBITS or ACTIVATES each of its target genes. Produce:
1. A research doc per drug at `docs/directionality/<DRUG_NAME>.md`
2. A single annotation JSON at `data/directionality_annotations.json`

## Candidates and their target genes

| Drug | Genes to research |
|---|---|
| LENALIDOMIDE | TNF, VEGFA, KDR, CTNNB1 |
| RITUXIMAB | IL6, JAK2, TGFB1 |
| IMATINIB | EGFR, JAK2, SMAD3, CTNNB1 |
| CELECOXIB | STAT3, PIK3CA, VEGFA, CTNNB1 |
| PEGINTERFERON ALFA-2B | EGFR, IL6, JAK2 |
| DEXAMETHASONE | SMAD3, VEGFA, KDR, CTNNB1 |
| THALIDOMIDE | TNF, VEGFA, CTNNB1 |
| SIROLIMUS | EGFR, JAK2, MTOR, PIK3CA |
| ASPIRIN | PIK3CA, TGFB1, CYP24A1 |
| CUCURBITACIN B | STAT3, SMAD3 |
| MILTEFOSINE | STAT3, SMAD3 |
| TRIFLUOPERAZINE | STAT3, SMAD3 |

## Background context

21 of 22 target genes are **pro_keloid** — they drive keloid formation through overexpression or overactivity. The exception is **TNF**, which is **context_dependent**: pro-proliferative via NF-κB but also induces MMPs that can degrade collagen. In keloid tissue, the anti-fibrotic MMP effect is largely neutralized by TIMP-1 overexpression — but the classification is context_dependent, so TNF drug-gene pairs should be marked `net_effect: "unclear"` regardless of inhibitor/activator status.

For pro_keloid genes, the question is: does the drug **inhibit** this gene's activity (beneficial) or **activate** it (unfavorable)?

## Research protocol per drug

For each drug, follow this sequence:

### Step 1: Query ChEMBL REST API

Fetch the drug's ChEMBL mechanisms. This gives structured action_type data.

1. Resolve drug name to ChEMBL ID:
   `GET https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name=<DRUG_NAME>&format=json`
   Extract `molecule_chembl_id` from `molecules[0]`.

2. Fetch mechanisms:
   `GET https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id=<ID>&format=json&limit=100`
   For each mechanism, note: `action_type`, `target_name`, `mechanism_of_action`.

3. Map ChEMBL target names to our gene symbols by matching on gene name or target description. `action_type` values: INHIBITOR → "inhibitor", ACTIVATOR → "activator", MODULATOR → "modulator", others → "unknown".

**Note on biologic candidates:** ChEMBL mechanisms are listed against the drug's primary known target. For biologics (RITUXIMAB, PEGINTERFERON ALFA-2B), the primary target (e.g., CD20/MS4A1 for rituximab) may not appear in our 22-gene list. In these cases, ChEMBL will return no matching gene — proceed directly to Step 2 (web fallback). This is expected and is not an error.

### Step 2: Web fallback for unknown pairs

For any drug-gene pair where ChEMBL returns no match or "unknown":

1. Browse DrugBank: search `https://go.drugbank.com/unearth/q?query=<drug_name>` and find the drug's mechanism of action page. Look for language about whether the drug inhibits or activates the target gene.

2. If still unclear, search PubMed: `https://pubmed.ncbi.nlm.nih.gov/?term=<drug_name>+<gene_symbol>+mechanism`. Read the top 3 abstracts for mechanism language.

3. If still unclear after both sources, mark as "unknown" with confidence "low".

### Step 3: Write research doc

Save to `docs/directionality/<DRUG_NAME>.md`:

```markdown
# <DRUG_NAME> — Directionality Research

**Drug class:** [brief description]
**ChEMBL ID:** [ID or "not found"]
**Research date:** [today]

## Gene-by-gene direction

| Gene | Pathway | Action | Net effect | Confidence | Source |
|---|---|---|---|---|---|
| MTOR | PI3K/AKT/mTOR | inhibitor | beneficial | high | ChEMBL:CHEMBL1614345 |
| ... | ... | ... | ... | ... | ... |

## Summary

[2-3 sentences: overall directionality verdict for this drug in keloid context. Is the drug's mechanism net-beneficial, net-unfavorable, or mixed?]

## Key sources

- [ChEMBL mechanism page or API result]
- [DrugBank URL if used]
- [PubMed PMID if used]
```

### Step 4: Append to annotation JSON

After each drug, append its entries to `data/directionality_annotations.json`. At the end of all 12 drugs, write the final JSON array. Schema per entry:

```json
{
  "drug_name": "SIROLIMUS",
  "gene_symbol": "MTOR",
  "action_direction": "inhibitor",
  "confidence": "high",
  "source": "ChEMBL:CHEMBL1614345",
  "keloid_target_role": "pro_keloid",
  "net_effect": "beneficial",
  "notes": "Sirolimus (rapamycin) is a canonical mTOR inhibitor. mTOR inhibition reduces fibroblast proliferation in keloid tissue."
}
```

`action_direction`: one of `inhibitor`, `activator`, `modulator`, `unknown`
`net_effect`: `beneficial` (inhibitor of pro_keloid gene), `unfavorable` (activator of pro_keloid gene), `unclear`
`confidence`: `high` (structured DB data), `medium` (clear web source), `low` (inferred or ambiguous)

## Output checklist

Before finishing, verify:
- [ ] 12 research docs created in `docs/directionality/`
- [ ] `data/directionality_annotations.json` is valid JSON array
- [ ] Every drug-gene pair has an entry (even if action_direction is "unknown")
- [ ] Every entry has all 7 fields
