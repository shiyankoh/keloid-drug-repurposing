# Keloid Pipeline — Deferred Work

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

## Every Cure MeDIC Database Check

**What:** Search Every Cure's open-source MeDIC database for keloid or fibrosis entries.

**Why:** Every Cure (Fajgenbaum's org) may have already scored drug-keloid repurposing candidates. Could validate pipeline output or surface candidates we missed.

**Effort:** 5 minutes, manual web search. No code needed.
