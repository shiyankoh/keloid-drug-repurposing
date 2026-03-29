import pytest
from scripts.db import init_db
from scripts.score_candidates import (
    query_overlaps,
    compute_scores,
    normalize_overlap_counts,
)


@pytest.fixture
def populated_db(db_conn):
    """DB with keloid targets and drugs that have known overlaps."""
    init_db(db_conn)

    # Insert keloid targets across 3 pathways
    db_conn.executemany(
        """INSERT INTO keloid_targets (gene_symbol, target_name, pathway, evidence_type, evidence_strength, source)
           VALUES (?, ?, ?, ?, ?, ?)""",
        [
            ("MTOR", "mTOR kinase", "PI3K/AKT/mTOR", "functional", 0.7, "manual"),
            ("PIK3CA", "PI3K alpha", "PI3K/AKT/mTOR", "functional", 0.6, "manual"),
            ("EGFR", "EGF receptor", "EGFR", "expression", 0.8, "manual"),
            ("TGFB1", "TGF-beta 1", "TGF-β/Smad", "functional", 0.9, "manual"),
            ("JAK2", "Janus kinase 2", "JAK/STAT", "functional", 0.5, "opentargets"),
        ],
    )

    # Drug A: hits 3 pathways (MTOR, EGFR, TGFB1) — should get multi-target bonus
    db_conn.execute(
        "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (1, 'DRUG_A', 'drug_a', 'approved')"
    )
    db_conn.executemany(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (?, ?, ?, ?)",
        [
            (1, "MTOR", "inhibitor", "DGIdb"),
            (1, "EGFR", "inhibitor", "DGIdb"),
            (1, "TGFB1", "inhibitor", "DGIdb"),
        ],
    )

    # Drug B: hits 1 pathway (MTOR only) — no multi-target bonus
    db_conn.execute(
        "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (2, 'DRUG_B', 'drug_b', 'approved')"
    )
    db_conn.execute(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (2, 'MTOR', 'inhibitor', 'DGIdb')"
    )

    # Drug C: hits 2 genes in SAME pathway (MTOR, PIK3CA both in PI3K/AKT/mTOR) — 1 distinct pathway
    db_conn.execute(
        "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (3, 'DRUG_C', 'drug_c', 'approved')"
    )
    db_conn.executemany(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (?, ?, ?, ?)",
        [
            (3, "MTOR", "inhibitor", "DGIdb"),
            (3, "PIK3CA", "inhibitor", "DGIdb"),
        ],
    )

    # Drug D: targets JAK2 only (low evidence pathway)
    db_conn.execute(
        "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (4, 'DRUG_D', 'drug_d', 'approved')"
    )
    db_conn.execute(
        "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (4, 'JAK2', 'inhibitor', 'DGIdb')"
    )

    db_conn.commit()
    return db_conn


class TestQueryOverlaps:
    def test_returns_all_matching_drugs(self, populated_db):
        overlaps = query_overlaps(populated_db)
        drug_names = set(o["drug_name"] for o in overlaps)
        assert drug_names == {"DRUG_A", "DRUG_B", "DRUG_C", "DRUG_D"}

    def test_drug_a_has_three_overlaps(self, populated_db):
        overlaps = query_overlaps(populated_db)
        drug_a = [o for o in overlaps if o["drug_name"] == "DRUG_A"]
        assert len(drug_a) == 3

    def test_no_overlaps_when_gene_not_in_keloid_targets(self, populated_db):
        populated_db.execute(
            "INSERT INTO drugs (id, drug_name, generic_name, approval_status) VALUES (99, 'ORPHAN', 'orphan', 'approved')"
        )
        populated_db.execute(
            "INSERT INTO drug_targets (drug_id, gene_symbol, action_type, source) VALUES (99, 'NONEXISTENT', 'inhibitor', 'DGIdb')"
        )
        populated_db.commit()
        overlaps = query_overlaps(populated_db)
        orphan = [o for o in overlaps if o["drug_name"] == "ORPHAN"]
        assert orphan == []


class TestComputeScores:
    def test_drug_a_gets_multi_target_bonus(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_a = [s for s in scores if s["drug_name"] == "DRUG_A"][0]
        assert drug_a["multi_target_bonus"] == 1.0

    def test_drug_b_no_multi_target_bonus(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_b = [s for s in scores if s["drug_name"] == "DRUG_B"][0]
        assert drug_b["multi_target_bonus"] == 0.0

    def test_drug_c_counts_one_distinct_pathway(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_c = [s for s in scores if s["drug_name"] == "DRUG_C"][0]
        assert drug_c["pathway_overlap_count"] == 1

    def test_drug_a_ranked_above_drug_b(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        scores.sort(key=lambda s: s["composite_score"], reverse=True)
        names = [s["drug_name"] for s in scores]
        assert names.index("DRUG_A") < names.index("DRUG_B")

    def test_avg_evidence_strength(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_b = [s for s in scores if s["drug_name"] == "DRUG_B"][0]
        assert drug_b["avg_evidence_strength"] == 0.7

    def test_drug_d_low_evidence(self, populated_db):
        overlaps = query_overlaps(populated_db)
        scores = compute_scores(overlaps)
        drug_d = [s for s in scores if s["drug_name"] == "DRUG_D"][0]
        assert drug_d["avg_evidence_strength"] == 0.5


class TestNormalizeOverlapCounts:
    def test_max_overlap_normalized_to_one(self):
        scores = [
            {"drug_name": "A", "pathway_overlap_count": 4},
            {"drug_name": "B", "pathway_overlap_count": 2},
            {"drug_name": "C", "pathway_overlap_count": 1},
        ]
        normalized = normalize_overlap_counts(scores)
        assert normalized[0]["normalized_overlap"] == 1.0
        assert normalized[1]["normalized_overlap"] == 0.5
        assert normalized[2]["normalized_overlap"] == 0.25

    def test_single_drug_normalized_to_one(self):
        scores = [{"drug_name": "A", "pathway_overlap_count": 1}]
        normalized = normalize_overlap_counts(scores)
        assert normalized[0]["normalized_overlap"] == 1.0

    def test_empty_list(self):
        assert normalize_overlap_counts([]) == []
