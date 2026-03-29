# Keloid Pipeline — Deferred Work

## Phase 2: Literature Corpus + Semantic Search

**What:** Ingest ~5K PubMed abstracts on keloids, embed with OpenAI `text-embedding-3-small`, store in Supabase with pgvector. Enable semantic search for drug-keloid hypotheses not captured by structured databases.

**Why:** Structured databases (OpenTargets, DGIdb) miss case reports and emerging research. Semantic search catches things like "sirolimus reduced keloid recurrence" in a case report abstract that has no structured drug-gene annotation.

**Pros:** Catches signal that structured DBs miss. Enables hypothesis-driven queries ("find papers about mTOR inhibitors and fibroblast proliferation").

**Cons:** Requires migration from SQLite to Supabase (pgvector). OpenAI API costs ~$5-10 for 5K embeddings. Significant infra work (embedding pipeline, pgvector index tuning).

**Depends on:** v1 complete. Migration to Supabase.

---

## Phase 2: Clinical Trials Cross-Reference

**What:** Query ClinicalTrials.gov API (v2, open, no auth) for any trials mentioning keloid/hypertrophic scar + drugs from the v1 candidate list. Store in a `clinical_trials` table.

**Why:** Validates candidates. If sirolimus already has a Phase 2 keloid trial, that's strong signal. Also identifies gaps — strong mechanistic rationale but zero trials = opportunity worth flagging to a dermatologist.

**Pros:** Direct validation signal. ClinicalTrials.gov API is simple and open. High value per line of code.

**Cons:** One additional script, small effort.

**Depends on:** v1 complete (need the candidate list to know what to search for).

---

## Phase 2: FDA FAERS Safety Profile Layer

**What:** Pull adverse event data from FDA FAERS for top candidates. Flag drugs with serious safety concerns that would make off-label use impractical.

**Why:** A drug that hits 5 keloid pathways but causes organ damage isn't a real candidate for off-label use. This is a practical filter for any drug you'd actually discuss with a dermatologist.

**Pros:** Critical for any candidate that moves beyond "interesting hypothesis" to "worth discussing with a doctor."

**Cons:** FAERS data is messy (inconsistent drug naming, duplicate reports). Requires non-trivial cleanup. Moderate effort.

**Depends on:** v1 complete.

---

## Pre-Pipeline: Every Cure MeDIC Database Check

**What:** Search Every Cure's open-source MeDIC (Mechanistic Disease-drug Interaction Candidates) database for keloid or fibrosis entries.

**Why:** Every Cure (Fajgenbaum's org) may have already scored drug-keloid repurposing candidates using their full bioinformatics pipeline. If so, you can validate your v1 output against their results — or skip building parts of the pipeline entirely.

**Pros:** Immediate signal. Zero code needed (web search / manual lookup). Should be done BEFORE building the pipeline.

**Cons:** None.

**Depends on:** Nothing. Do this first.
