# Keloid Pipeline — Deferred Work

## Science Deepening Roadmap (priority order)

Four gaps identified 2026-04-15. Tackling in this order:

- [ ] **A. Directionality** — can't tell inhibitors from activators. Current status: IN PROGRESS (Tasks 0–6 done. Task 7 RUNNING overnight 2026-04-15 via `ralph/directionality-research/ralph.sh`.)
- [ ] **B. Evidence quality** — most evidence scores are a flat 0.70. Scoring is too uniform to distinguish strongly vs weakly supported candidates.
- [ ] **C. Missing data sources** — ChEMBL (potency/IC50), GEO (actual keloid gene expression), DrugBank (mechanism detail) not yet tapped.
- [ ] **D. Literature blind spot** — structured DBs miss case reports and emerging research. Need PubMed semantic search.

---

## Directionality Layer — Task 7 RUNNING

**Status (2026-04-15 evening):** Ralph loop launched in a separate terminal against `ralph/directionality-research/prd.json` (14 stories: init + 12 drugs + final validation). ETA 1-2 hours.

**Files:**
- `docs/ralph-prd-directionality.md` — human-readable PRD (committed)
- `ralph/directionality-research/prd.json` — 14 user stories for ralph (committed)
- `ralph/directionality-research/ralph.sh` + `CLAUDE.md` — copied from ralph-skills plugin (gitignored)
- Ralph runs on branch `ralph/directionality-research` (auto-created from main)

**Monitor:**
```bash
tail -f ralph/directionality-research/progress.txt
ls docs/directionality/   # grows to 12 .md files
```

**When ralph completes (`<promise>COMPLETE</promise>`):**
1. Review 12 research docs in `docs/directionality/`
2. Spot-check `data/directionality_annotations.json` (expect 38 entries, 8 fields each)
3. Merge ralph branch into main: `git checkout main && git merge ralph/directionality-research`
4. Import + rescore + regenerate report:
```bash
source venv/bin/activate
python3 -m scripts.08_import_directionality
python3 -m scripts.03_score_candidates
python3 -m scripts.04_generate_report
```
5. Report will show `[+]/[!]/[~]` instead of `[?]` for the 12 candidates

**If ralph gets stuck / fails a story:** Edit `prd.json` to set that story's `passes: true` manually if the output is good enough, or adjust acceptance criteria and rerun `./ralph.sh --tool claude N`.

---

## DONE: Clinical Trials Cross-Reference

Completed 2026-03-29. See `scripts/clinical_trials.py`. Found 269 keloid/scar trials, matched 6 to top 100 candidates. Results integrated into report.

---

## DONE: Severity Tiers

Completed 2026-03-29. See `data/severity_tiers.csv`. Top 50 drugs classified as green/yellow/red. Report leads with actionable candidates (green/yellow, 2+ pathways).

---

## Literature Corpus + Semantic Search

**What:** Ingest ~5K PubMed abstracts on keloids, embed with OpenAI `text-embedding-3-small`, store in Supabase with pgvector. Enable semantic search for drug-keloid hypotheses not captured by structured databases.

**Why:** Structured databases (OpenTargets, DGIdb) miss case reports and emerging research. Semantic search catches things like "sirolimus reduced keloid recurrence" in a case report abstract that has no structured drug-gene annotation.

**Pros:** Catches signal that structured DBs miss. Enables hypothesis-driven queries ("find papers about mTOR inhibitors and fibroblast proliferation").

**Cons:** Requires migration from SQLite to Supabase (pgvector). OpenAI API costs ~$5-10 for 5K embeddings. Significant infra work (embedding pipeline, pgvector index tuning).

**Depends on:** Migration to Supabase.

---

## FDA FAERS Safety Profile Layer

**What:** Pull adverse event data from FDA FAERS for top candidates. Flag drugs with serious safety concerns that would make off-label use impractical.

**Why:** Partially addressed by manual severity tiers, but FAERS data would add quantitative adverse event frequencies and catch drugs outside the top 50 that were manually classified.

**Cons:** FAERS data is messy (inconsistent drug naming, duplicate reports). Moderate effort. Lower priority now that severity tiers exist.

---

## Known Limitations (Clinical Review)

The following limitations were flagged by a physician reviewing the pipeline output. These aren't bugs — they're fundamental gaps in the computational approach that should inform how results are interpreted and what future work could address.

### 1. Molecular interaction directionality ignored

Inhibitors, activators, and partial agonists are treated equally. A drug that activates a keloid-promoting pathway is scored the same as one that inhibits it. Worse, some drugs have tissue-specific effects — e.g., tamoxifen suppresses estrogen in breast tissue but activates it in endometrial tissue. The pipeline has no way to distinguish "hits the pathway" from "hits the pathway in the right direction in skin fibroblasts."

### 2. Route of delivery and pharmacokinetics not modeled

Topical, intralesional, and systemic delivery have very different efficacy and side-effect profiles. Intralesional steroids work well for keloids; systemic steroids would be inappropriate. The pipeline treats all drugs as if route doesn't matter, which overstates the viability of some candidates and understates others.

### 3. Effect size and pathway importance not weighted

All pathways are treated equally, but clinically some targets matter far more than others for keloids. Keloids are likely driven primarily by pro-inflammatory processes (which is why intralesional steroids work and why so many anti-inflammatory drugs score highly). mTOR, for example, is strongly implicated in vascular malformations but probably plays a minor role in keloid pathogenesis. The pipeline rewards breadth of pathway coverage rather than depth on the pathways that actually drive keloid biology.

### 4. Mechanobiology ignored

From a surgical perspective, mechanical tension and scar placement are the biggest factors in keloid prevention and the basis of scar revision. Silicone taping works because it modifies tension, not because of a molecular pathway. Drug delivery mechanisms (microneedling, intralesional depots) that exploit mechanical properties aren't captured. This is a fundamentally different treatment axis that a molecular pipeline can't model.

### 5. Wound healing temporal dynamics not captured

Keloid formation is time-dependent and occurs in discrete phases. Different therapies are indicated at different phases (acute inflammation, proliferation, remodeling). The pipeline treats keloid formation as a static state rather than a process, so it can't distinguish drugs that should be used early (to prevent keloid formation) from those that address established keloids.

---

## Every Cure MeDIC Database Check

**What:** Search Every Cure's open-source MeDIC database for keloid or fibrosis entries.

**Why:** Every Cure (Fajgenbaum's org) may have already scored drug-keloid repurposing candidates. Could validate pipeline output or surface candidates we missed.

**Effort:** 5 minutes, manual web search. No code needed.
