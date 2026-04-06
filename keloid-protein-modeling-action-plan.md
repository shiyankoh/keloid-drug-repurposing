# Protein-Substrate Modeling for Keloid Drug Candidate Validation

An action plan to add structural and computational validation to drug repurposing candidates, confirming that a candidate drug can physically bind to and modulate its keloid-relevant target protein.

---

## Why This Matters

Pathway overlap scoring tells you *a drug is known to affect a pathway involved in keloids*. Protein-substrate modeling tells you *how well the drug molecule actually fits into the target protein's binding site*. This is the difference between "this key might fit that lock" and "here's a 3D simulation showing the key turning."

This layer filters out false positives where a drug nominally targets the right pathway but has weak or implausible binding to the specific protein isoform relevant in keloid fibroblasts.

---

## What We're Modeling

For each top drug candidate from the pipeline:

1. **Retrieve the 3D structure** of the keloid-relevant target protein
2. **Retrieve or generate the 3D structure** of the drug molecule
3. **Dock the drug into the protein's binding site** and score the interaction
4. **Estimate binding affinity** — does it bind tightly enough to have a therapeutic effect?
5. **Compare against known binders** — how does the candidate compare to drugs already proven to modulate this target?

---

## Phase 1: Structural Data Acquisition (Week 1)

### 1A. Protein Structures

**Primary source:** RCSB Protein Data Bank (PDB)
- URL: rcsb.org / API: data.rcsb.org
- Search for crystal/cryo-EM structures of each keloid target protein
- Priority targets with likely available structures:

| Target | Gene | Role in Keloids | PDB Likely? |
|--------|------|-----------------|-------------|
| TGF-β receptor | TGFBR1/TGFBR2 | Master fibrosis driver | Yes |
| mTOR | MTOR | Fibroblast proliferation | Yes |
| PI3K | PIK3CA | Proliferation/survival signaling | Yes |
| AKT | AKT1 | Downstream of PI3K | Yes |
| IL-6 receptor | IL6R | Inflammatory signaling | Yes |
| EGFR | EGFR | Growth signaling | Yes |
| JAK1/2 | JAK1/JAK2 | Cytokine signaling | Yes |
| Prolyl-tRNA synthetase | PARS2/EPRS1 | Collagen biosynthesis | Partial |
| CYP24A1 | CYP24A1 | Vitamin D metabolism | Yes |
| COL1A1/COL1A2 | COL1A1 | Collagen itself | Structural only |

**Fallback for missing structures:** AlphaFold Protein Structure Database
- URL: alphafold.ebi.ac.uk
- Covers nearly all human proteins with high-confidence predicted structures
- API: alphafold.ebi.ac.uk/api — query by UniProt ID

**Action items:**
- [ ] For each target in `keloid_targets` table, query PDB API for available structures
- [ ] Where no experimental structure exists, pull AlphaFold prediction
- [ ] Store in a `protein_structures` table: `gene_symbol`, `pdb_id`, `alphafold_id`, `structure_source`, `resolution_angstroms`, `file_path`
- [ ] Download .pdb or .cif files to local storage

### 1B. Drug (Ligand) Structures

**Primary sources:**
- **PubChem:** pubchem.ncbi.nlm.nih.gov — SDF/MOL files for all approved drugs
- **DrugBank:** 3D SDF structures included in data downloads
- **ZINC20:** zinc20.docking.org — pre-prepared ligand files ready for docking

**Action items:**
- [ ] For each drug in `drugs` table, pull 3D structure (SDF format) from PubChem API
- [ ] Convert to PDBQT format (required by most docking software) using Open Babel
- [ ] Store in `drug_structures` table: `drug_id`, `pubchem_cid`, `smiles`, `file_path`

---

## Phase 2: Molecular Docking Setup (Week 2)

### Tool Selection

| Tool | Type | Difficulty | Best For |
|------|------|-----------|----------|
| **AutoDock Vina** | Classical docking | Low | Fast screening of many candidates |
| **DiffDock** | ML-based docking | Medium | Higher accuracy, no binding site needed |
| **Uni-Mol** | ML-based | Medium | Binding affinity prediction |
| **P2Rank** | Binding site predictor | Low | Finding where drug should bind |
| **PLIP** | Interaction profiler | Low | Analyzing what holds the drug in place |

**Recommended primary stack:** AutoDock Vina for speed + DiffDock for accuracy on top candidates.

### 2A. AutoDock Vina Pipeline (Bulk Screening)

AutoDock Vina is the workhorse of computational drug screening — fast, well-documented, and scriptable.

```bash
# Install
pip install vina
# or
conda install -c conda-forge autodock-vina

# Also need:
pip install meeko          # Ligand preparation
pip install openbabel      # Format conversion
```

**Pipeline per candidate pair (protein + drug):**

1. **Prepare protein:** Remove water molecules, add hydrogens, convert to PDBQT
2. **Identify binding site:** Use known co-crystallized ligand position, or run P2Rank
3. **Prepare ligand:** Convert drug SDF → PDBQT, generate conformers
4. **Define search box:** Centered on binding site, typically 20–30 Å cube
5. **Run Vina:** Returns binding poses ranked by predicted binding energy (kcal/mol)
6. **Interpret:** More negative = stronger binding. Threshold: < -7.0 kcal/mol is promising

```python
# Pseudocode for batch docking
from vina import Vina

def dock_candidate(protein_pdbqt, ligand_pdbqt, center, box_size):
    v = Vina(sf_name='vina')
    v.set_receptor(protein_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=32, n_poses=10)
    return v.energies(n_poses=5)  # Top 5 poses with scores
```

### 2B. DiffDock (ML-Based, Higher Accuracy)

DiffDock uses diffusion models to predict binding poses without needing a predefined binding site — useful when you're not sure exactly where the drug binds.

- **Repo:** github.com/gcorso/DiffDock
- **Runs on:** GPU preferred, CPU possible but slow
- **Input:** Protein .pdb + drug SMILES string
- **Output:** Predicted binding poses with confidence scores

Reserve DiffDock for your top 10–15 candidates from the Vina screen.

### 2C. Binding Site Prediction with P2Rank

When you don't have a co-crystallized ligand to define the binding pocket:

- **Tool:** P2Rank (github.com/rdk/p2rank)
- **Input:** Protein .pdb file
- **Output:** Ranked list of predicted binding pockets with coordinates
- **Use case:** Feed pocket coordinates into Vina's search box definition

---

## Phase 3: Binding Affinity & Interaction Analysis (Week 3)

### 3A. Binding Affinity Estimation

Docking scores from Vina are rough proxies for binding affinity. For more confidence:

**Option 1: MM-GBSA Rescoring**
- Post-process Vina poses with molecular mechanics / generalized Born surface area
- More physically accurate than docking scores alone
- Tool: AmberTools (free) or OpenMM

**Option 2: ML-Based Affinity Prediction**
- **Uni-Mol / OnionNet-2:** Feed docking poses → predicted Kd/Ki values
- Faster than MM-GBSA, increasingly competitive in accuracy
- Good for ranking candidates relative to each other

**Practical threshold guidance:**
- Kd < 100 nM: Strong binder, likely pharmacologically active at therapeutic doses
- Kd 100 nM – 1 μM: Moderate, could work at higher doses
- Kd > 10 μM: Weak, unlikely to be therapeutically relevant

### 3B. Interaction Profiling with PLIP

PLIP (Protein-Ligand Interaction Profiler) tells you *why* a drug binds:

- **URL:** plip-tool.biotec.tu-dresden.de (web) or install locally
- **Input:** Docked protein-ligand complex (.pdb)
- **Output:** Catalog of interactions — hydrogen bonds, hydrophobic contacts, salt bridges, π-stacking
- **Value:** A drug making 4+ hydrogen bonds with catalytic residues is more convincing than one held in place by a single hydrophobic contact

### 3C. Comparative Analysis

For each keloid target, compare your candidate's docking score against:
1. **Known active drugs** for that target (positive controls)
2. **Random approved drugs** not expected to bind (negative controls)
3. **Drugs already used off-label for keloids** (benchmark)

If your candidate scores comparably to known actives, that's strong validation.

---

## Phase 4: Integration with Main Pipeline (Week 4)

### Data Model Extension

```
protein_structures
├── id (uuid)
├── gene_symbol (text)
├── pdb_id (text, nullable)
├── alphafold_id (text, nullable)
├── structure_source (text)  -- 'pdb_experimental' | 'alphafold_predicted'
├── resolution_angstroms (float, nullable)
├── binding_site_coords (jsonb)  -- from P2Rank or co-crystal
└── file_path (text)

drug_structures
├── id (uuid)
├── drug_id (uuid → drugs)
├── pubchem_cid (text)
├── smiles (text)
├── file_path (text)
└── format (text)  -- 'sdf' | 'pdbqt' | 'mol2'

docking_results
├── id (uuid)
├── drug_id (uuid → drugs)
├── protein_id (uuid → protein_structures)
├── method (text)  -- 'vina' | 'diffdock'
├── binding_energy_kcal (float)
├── confidence_score (float)
├── predicted_kd_nm (float, nullable)
├── interaction_count (int)  -- from PLIP
├── interaction_types (jsonb)  -- {"h_bonds": 3, "hydrophobic": 5, ...}
├── pose_file_path (text)
├── vs_known_active_percentile (float)  -- how it ranks vs positive controls
└── created_at (timestamp)
```

### Updated Composite Score

Revise the main pipeline scoring to incorporate structural evidence:

```
candidate_score = (
    pathway_overlap_count      × 0.20 +
    evidence_strength_avg      × 0.20 +
    structural_binding_score   × 0.20 +  ← NEW
    literature_mention_score   × 0.15 +
    clinical_trial_exists      × 0.10 +
    safety_profile_score       × 0.10 +
    multi_target_bonus         × 0.05
)
```

Where `structural_binding_score` is normalized from docking energy, scaled by:
- Comparison to known active compounds for that target
- Number and quality of protein-ligand interactions
- Confidence of the protein structure itself (experimental > AlphaFold)

---

## Compute Requirements

| Step | CPU Only | With GPU | Cloud Option |
|------|----------|----------|--------------|
| AutoDock Vina (per pair) | 2–10 min | N/A (CPU tool) | — |
| Vina batch (100 pairs) | 3–15 hrs | N/A | AWS Batch |
| DiffDock (per pair) | 20–60 min | 2–5 min | Colab/Lambda |
| DiffDock batch (15 pairs) | 5–15 hrs | 30–75 min | Colab/Lambda |
| P2Rank (per protein) | < 1 min | N/A | — |
| PLIP (per complex) | < 1 min | N/A | — |
| MM-GBSA (per complex) | 30–120 min | 5–15 min | — |

Your M5 MacBook Pro can handle Vina screening and P2Rank locally. DiffDock benefits from GPU — consider Google Colab Pro ($10/mo) or a Lambda Labs spot instance for the accuracy pass.

---

## Practical Workflow Summary

```
1. Pull protein structures (PDB / AlphaFold)
          ↓
2. Pull drug structures (PubChem)
          ↓
3. Predict binding sites (P2Rank)
          ↓
4. Bulk screen with AutoDock Vina (~100 drug-protein pairs)
          ↓
5. Filter: keep candidates with binding energy < -7.0 kcal/mol
          ↓
6. Re-dock top 15 with DiffDock for higher confidence
          ↓
7. Profile interactions with PLIP
          ↓
8. Estimate binding affinity (Uni-Mol or MM-GBSA)
          ↓
9. Compare against known actives (percentile ranking)
          ↓
10. Feed scores back into main pipeline composite score
```

---

## Key Risks & Mitigations

| Risk | Mitigation |
|------|-----------|
| Docking scores don't always predict real-world activity | Use as a ranking/filtering tool, not a binary yes/no. Validate top candidates against literature case reports. |
| AlphaFold structures may have low confidence in binding site regions | Check per-residue confidence (pLDDT). Only use regions with pLDDT > 70. Flag structures below threshold. |
| Drug may bind the right protein but in the wrong tissue context | Cross-reference with keloid fibroblast gene expression data (GEO database) to confirm the target is actually expressed in keloid tissue. |
| Conformational dynamics matter — static docking misses flexibility | For top 3–5 candidates, consider running short molecular dynamics simulations (OpenMM) to assess binding stability. This is a stretch goal. |
| You are not a computational chemist | This pipeline produces ranked hypotheses to discuss with a dermatologist and potentially a computational biology collaborator. It does not produce prescriptions. |

---

## Resources to Get Started

- **AutoDock Vina docs:** autodock-vina.readthedocs.io
- **DiffDock paper:** Corso et al., ICLR 2023 — "DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking"
- **P2Rank:** github.com/rdk/p2rank
- **PLIP web:** plip-tool.biotec.tu-dresden.de
- **AlphaFold DB API:** alphafold.ebi.ac.uk/api-docs
- **PubChem API:** pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- **Open Babel:** openbabel.org (format conversion Swiss army knife)
- **Practical tutorial:** "Molecular Docking with AutoDock Vina" on TeachOpenCADD (github.com/volkamerlab/teachopencadd)
