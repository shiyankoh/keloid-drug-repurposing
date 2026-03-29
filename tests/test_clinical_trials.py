import json
import os
import pytest
from scripts.clinical_trials import (
    parse_trials_response,
    match_trials_to_candidates,
)


@pytest.fixture
def trials_response():
    fixture_path = os.path.join(
        os.path.dirname(__file__), "fixtures", "clinicaltrials_keloid.json"
    )
    with open(fixture_path) as f:
        return json.load(f)


@pytest.fixture
def candidate_names():
    return ["SIROLIMUS", "IMATINIB", "VERAPAMIL", "CELECOXIB", "DEXAMETHASONE"]


class TestParseTrialsResponse:
    def test_extracts_all_trials(self, trials_response):
        trials = parse_trials_response(trials_response)
        assert len(trials) == 4

    def test_extracts_nct_id(self, trials_response):
        trials = parse_trials_response(trials_response)
        nct_ids = [t["nct_id"] for t in trials]
        assert "NCT04049552" in nct_ids

    def test_extracts_phase(self, trials_response):
        trials = parse_trials_response(trials_response)
        rapa = [t for t in trials if t["nct_id"] == "NCT04049552"][0]
        assert rapa["phase"] == "EARLY_PHASE1"

    def test_extracts_status(self, trials_response):
        trials = parse_trials_response(trials_response)
        rapa = [t for t in trials if t["nct_id"] == "NCT04049552"][0]
        assert rapa["status"] == "COMPLETED"

    def test_extracts_intervention_names(self, trials_response):
        trials = parse_trials_response(trials_response)
        rapa = [t for t in trials if t["nct_id"] == "NCT04049552"][0]
        assert "Rapamycin 8% Ointment" in rapa["interventions"]
        assert "Placebo" in rapa["interventions"]

    def test_extracts_conditions(self, trials_response):
        trials = parse_trials_response(trials_response)
        imatinib = [t for t in trials if t["nct_id"] == "NCT00000001"][0]
        assert "Keloid" in imatinib["conditions"]
        assert "Hypertrophic Scar" in imatinib["conditions"]

    def test_empty_response(self):
        empty = {"studies": []}
        trials = parse_trials_response(empty)
        assert trials == []


class TestMatchTrialsToCandidates:
    def test_matches_sirolimus_via_rapamycin(self, trials_response, candidate_names):
        """Sirolimus is also known as rapamycin. Should match the RAPA-Keloid trial."""
        trials = parse_trials_response(trials_response)
        matches = match_trials_to_candidates(trials, candidate_names)
        assert "SIROLIMUS" in matches
        assert any(t["nct_id"] == "NCT04049552" for t in matches["SIROLIMUS"])

    def test_matches_imatinib(self, trials_response, candidate_names):
        trials = parse_trials_response(trials_response)
        matches = match_trials_to_candidates(trials, candidate_names)
        assert "IMATINIB" in matches
        assert len(matches["IMATINIB"]) == 1

    def test_matches_verapamil(self, trials_response, candidate_names):
        trials = parse_trials_response(trials_response)
        matches = match_trials_to_candidates(trials, candidate_names)
        assert "VERAPAMIL" in matches

    def test_unmatched_candidate_not_in_results(self, trials_response, candidate_names):
        """CELECOXIB and DEXAMETHASONE have no trials in the fixture."""
        trials = parse_trials_response(trials_response)
        matches = match_trials_to_candidates(trials, candidate_names)
        assert "CELECOXIB" not in matches
        assert "DEXAMETHASONE" not in matches

    def test_mystery_drug_not_matched(self, trials_response, candidate_names):
        """MYSTERY COMPOUND XYZ doesn't match any candidate."""
        trials = parse_trials_response(trials_response)
        matches = match_trials_to_candidates(trials, candidate_names)
        matched_ncts = []
        for drug_trials in matches.values():
            matched_ncts.extend(t["nct_id"] for t in drug_trials)
        assert "NCT00000003" not in matched_ncts

    def test_empty_candidates(self, trials_response):
        trials = parse_trials_response(trials_response)
        matches = match_trials_to_candidates(trials, [])
        assert matches == {}

    def test_empty_trials(self, candidate_names):
        matches = match_trials_to_candidates([], candidate_names)
        assert matches == {}
