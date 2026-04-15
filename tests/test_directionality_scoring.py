import pytest
from scripts.score_candidates import compute_directionality_flag


def test_all_beneficial_is_favorable():
    gene_details = [
        {"gene_symbol": "MTOR", "action_direction": "inhibitor", "keloid_target_role": "pro_keloid"},
        {"gene_symbol": "EGFR", "action_direction": "inhibitor", "keloid_target_role": "pro_keloid"},
    ]
    assert compute_directionality_flag(gene_details) == "favorable"


def test_all_unfavorable_is_unfavorable():
    gene_details = [
        {"gene_symbol": "MTOR", "action_direction": "activator", "keloid_target_role": "pro_keloid"},
    ]
    assert compute_directionality_flag(gene_details) == "unfavorable"


def test_mixed_is_mixed():
    gene_details = [
        {"gene_symbol": "MTOR", "action_direction": "inhibitor", "keloid_target_role": "pro_keloid"},
        {"gene_symbol": "EGFR", "action_direction": "activator", "keloid_target_role": "pro_keloid"},
    ]
    assert compute_directionality_flag(gene_details) == "mixed"


def test_all_unknown_is_unknown():
    gene_details = [
        {"gene_symbol": "MTOR", "action_direction": "unknown", "keloid_target_role": "pro_keloid"},
        {"gene_symbol": "EGFR", "action_direction": None, "keloid_target_role": "pro_keloid"},
    ]
    assert compute_directionality_flag(gene_details) == "unknown"


def test_empty_gene_details_is_unknown():
    assert compute_directionality_flag([]) == "unknown"


def test_modulator_counts_as_unknown():
    gene_details = [
        {"gene_symbol": "MTOR", "action_direction": "modulator", "keloid_target_role": "pro_keloid"},
    ]
    assert compute_directionality_flag(gene_details) == "unknown"
