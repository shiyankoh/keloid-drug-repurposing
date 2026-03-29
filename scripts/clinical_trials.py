"""Cross-reference drug candidates against ClinicalTrials.gov for keloid/scar trials."""

import json
import logging
import os
import requests

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

CT_API_URL = "https://clinicaltrials.gov/api/v2/studies"
RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
SCORES_PATH = os.path.join(RAW_DIR, "scores.json")
OUTPUT_PATH = os.path.join(RAW_DIR, "clinical_trials.json")
PAGE_SIZE = 100

# Common drug name aliases (generic ↔ brand/chemical name)
# Used when matching trial intervention names to our candidate list
DRUG_ALIASES = {
    "SIROLIMUS": ["rapamycin", "sirolimus"],
    "FLUOROURACIL": ["5-fu", "fluorouracil", "5-fluorouracil"],
    "DEXAMETHASONE": ["dexamethasone", "decadron"],
    "VERAPAMIL": ["verapamil"],
    "BLEOMYCIN": ["bleomycin"],
    "IMATINIB": ["imatinib", "gleevec", "glivec"],
    "TAMOXIFEN": ["tamoxifen", "nolvadex"],
    "CELECOXIB": ["celecoxib", "celebrex"],
    "ASPIRIN": ["aspirin", "acetylsalicylic"],
    "INTERFERON": ["interferon"],
    "RITUXIMAB": ["rituximab", "rituxan", "mabthera"],
    "THALIDOMIDE": ["thalidomide"],
    "LENALIDOMIDE": ["lenalidomide", "revlimid"],
    "PIRFENIDONE": ["pirfenidone", "esbriet"],
    "NINTEDANIB": ["nintedanib", "ofev"],
    "COLCHICINE": ["colchicine"],
    "METHOTREXATE": ["methotrexate"],
}


def fetch_keloid_trials(cache_path=None):
    """Fetch all keloid and hypertrophic scar trials from ClinicalTrials.gov.
    Uses cache if available."""
    cache = cache_path or os.path.join(RAW_DIR, "clinicaltrials_keloid.json")

    if os.path.exists(cache):
        logger.info(f"Loading cached ClinicalTrials.gov response from {cache}")
        with open(cache) as f:
            return json.load(f)

    all_studies = []
    fields = "NCTId,BriefTitle,OverallStatus,Phase,StartDate,Condition,InterventionName"

    for condition in ["keloid", "hypertrophic scar"]:
        page_token = None
        while True:
            params = {
                "query.cond": condition,
                "pageSize": PAGE_SIZE,
                "fields": fields,
                "countTotal": "true",
            }
            if page_token:
                params["pageToken"] = page_token

            logger.info(f"Fetching trials for '{condition}'...")
            resp = requests.get(CT_API_URL, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()

            studies = data.get("studies", [])
            all_studies.extend(studies)

            page_token = data.get("nextPageToken")
            if not page_token:
                break

    # Deduplicate by NCT ID (a trial may appear under both conditions)
    seen = set()
    unique_studies = []
    for study in all_studies:
        nct = study["protocolSection"]["identificationModule"]["nctId"]
        if nct not in seen:
            seen.add(nct)
            unique_studies.append(study)

    result = {"studies": unique_studies, "totalCount": len(unique_studies)}

    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        json.dump(result, f, indent=2)
    logger.info(f"Cached {len(unique_studies)} trials to {cache}")

    return result


def parse_trials_response(response_json):
    """Parse ClinicalTrials.gov response into a flat list of trial dicts."""
    studies = response_json.get("studies", [])
    trials = []
    for study in studies:
        proto = study.get("protocolSection", {})
        ident = proto.get("identificationModule", {})
        status = proto.get("statusModule", {})
        conditions = proto.get("conditionsModule", {})
        design = proto.get("designModule", {})
        arms = proto.get("armsInterventionsModule", {})

        interventions = [i["name"] for i in arms.get("interventions", [])]
        phases = design.get("phases", [])

        trials.append({
            "nct_id": ident.get("nctId", ""),
            "title": ident.get("briefTitle", ""),
            "status": status.get("overallStatus", ""),
            "start_date": status.get("startDateStruct", {}).get("date", ""),
            "conditions": conditions.get("conditions", []),
            "phase": phases[0] if phases else "N/A",
            "interventions": interventions,
        })
    return trials


def match_trials_to_candidates(trials, candidate_names):
    """Match trials to candidate drug names by checking intervention names.
    Uses DRUG_ALIASES for common name variations.
    Returns dict of drug_name → [matched trials]."""
    if not trials or not candidate_names:
        return {}

    # Build lookup: alias_lowercase → canonical drug name
    alias_lookup = {}
    for drug_name in candidate_names:
        # Always add the drug name itself (lowercased)
        alias_lookup[drug_name.lower()] = drug_name
        # Add known aliases
        if drug_name in DRUG_ALIASES:
            for alias in DRUG_ALIASES[drug_name]:
                alias_lookup[alias.lower()] = drug_name

    matches = {}
    for trial in trials:
        intervention_text = " ".join(trial["interventions"]).lower()
        title_text = trial["title"].lower()
        search_text = intervention_text + " " + title_text

        for alias, canonical_name in alias_lookup.items():
            if alias in search_text:
                if canonical_name not in matches:
                    matches[canonical_name] = []
                # Avoid duplicate trial entries for same drug
                if not any(t["nct_id"] == trial["nct_id"] for t in matches[canonical_name]):
                    matches[canonical_name].append(trial)

    return matches


def main():
    # Load scored candidates
    if not os.path.exists(SCORES_PATH):
        logger.error(f"Scores not found at {SCORES_PATH}. Run 03_score_candidates.py first.")
        return

    with open(SCORES_PATH) as f:
        scores = json.load(f)

    # Use top 100 candidates for matching
    candidate_names = [s["drug_name"] for s in scores[:100]]

    # Fetch trials
    response = fetch_keloid_trials()
    trials = parse_trials_response(response)
    logger.info(f"Parsed {len(trials)} keloid/scar trials")

    # Match
    matches = match_trials_to_candidates(trials, candidate_names)
    logger.info(f"Found trial matches for {len(matches)} of {len(candidate_names)} candidates")

    # Add trial data back to scores
    for s in scores:
        drug_trials = matches.get(s["drug_name"], [])
        s["clinical_trials"] = drug_trials
        s["has_keloid_trial"] = len(drug_trials) > 0

    # Save updated scores
    with open(SCORES_PATH, "w") as f:
        json.dump(scores, f, indent=2, default=list)

    # Print summary
    print(f"\n{'Drug':<30}{'Trials':<8}{'Phases':<25}{'Status'}")
    print("-" * 80)
    for drug_name in sorted(matches.keys()):
        drug_trials = matches[drug_name]
        phases = ", ".join(sorted(set(t["phase"] for t in drug_trials)))
        statuses = ", ".join(sorted(set(t["status"] for t in drug_trials)))
        print(f"{drug_name:<30}{len(drug_trials):<8}{phases:<25}{statuses}")

    # Save trial matches separately for easy reference
    with open(OUTPUT_PATH, "w") as f:
        json.dump(matches, f, indent=2, default=list)
    logger.info(f"Trial matches saved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
