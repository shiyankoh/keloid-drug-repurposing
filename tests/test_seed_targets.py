import csv
import io
import json
import os
import pytest
from scripts.db import init_db
from scripts.seed_keloid_targets import (
    parse_seed_csv,
    parse_opentargets_response,
    insert_targets,
)


@pytest.fixture
def sample_csv():
    return """gene_symbol,target_name,pathway,evidence_type,evidence_strength,source
TGFB1,transforming growth factor beta 1,TGF-β/Smad,functional,0.7,manual_curation
MTOR,mechanistic target of rapamycin kinase,PI3K/AKT/mTOR,functional,0.7,manual_curation"""


@pytest.fixture
def opentargets_response():
    fixture_path = os.path.join(
        os.path.dirname(__file__), "fixtures", "opentargets_keloid.json"
    )
    with open(fixture_path) as f:
        return json.load(f)


class TestParseSeedCsv:
    def test_returns_list_of_dicts(self, sample_csv):
        targets = parse_seed_csv(io.StringIO(sample_csv))
        assert len(targets) == 2
        assert targets[0]["gene_symbol"] == "TGFB1"
        assert targets[1]["gene_symbol"] == "MTOR"

    def test_evidence_strength_is_float(self, sample_csv):
        targets = parse_seed_csv(io.StringIO(sample_csv))
        assert isinstance(targets[0]["evidence_strength"], float)
        assert targets[0]["evidence_strength"] == 0.7

    def test_empty_csv_returns_empty_list(self):
        empty = "gene_symbol,target_name,pathway,evidence_type,evidence_strength,source\n"
        targets = parse_seed_csv(io.StringIO(empty))
        assert targets == []


class TestParseOpentargetsResponse:
    def test_extracts_gene_symbols(self, opentargets_response):
        targets = parse_opentargets_response(opentargets_response)
        symbols = [t["gene_symbol"] for t in targets]
        assert "PHLDA3" in symbols
        assert "NEDD4" in symbols
        assert "TGFB1" in symbols

    def test_includes_ensembl_id(self, opentargets_response):
        targets = parse_opentargets_response(opentargets_response)
        phlda3 = [t for t in targets if t["gene_symbol"] == "PHLDA3"][0]
        assert phlda3["ensembl_id"] == "ENSG00000174307"

    def test_uses_opentargets_score(self, opentargets_response):
        targets = parse_opentargets_response(opentargets_response)
        phlda3 = [t for t in targets if t["gene_symbol"] == "PHLDA3"][0]
        assert phlda3["evidence_strength"] == 0.52

    def test_empty_response_returns_empty_list(self):
        empty = {"data": {"disease": None}}
        targets = parse_opentargets_response(empty)
        assert targets == []

    def test_zero_targets_returns_empty_list(self):
        zero = {"data": {"disease": {"associatedTargets": {"count": 0, "rows": []}}}}
        targets = parse_opentargets_response(zero)
        assert targets == []


class TestInsertTargets:
    def test_inserts_csv_targets(self, db_conn, sample_csv):
        init_db(db_conn)
        csv_targets = parse_seed_csv(io.StringIO(sample_csv))
        insert_targets(db_conn, csv_targets, [])
        cursor = db_conn.execute("SELECT COUNT(*) as cnt FROM keloid_targets")
        assert cursor.fetchone()["cnt"] == 2

    def test_csv_seeds_win_on_conflict(self, db_conn, sample_csv, opentargets_response):
        init_db(db_conn)
        csv_targets = parse_seed_csv(io.StringIO(sample_csv))
        api_targets = parse_opentargets_response(opentargets_response)
        insert_targets(db_conn, csv_targets, api_targets)
        cursor = db_conn.execute(
            "SELECT evidence_strength, source FROM keloid_targets WHERE gene_symbol = 'TGFB1'"
        )
        row = cursor.fetchone()
        assert row["evidence_strength"] == 0.7
        assert row["source"] == "manual_curation"

    def test_api_targets_added_when_no_conflict(self, db_conn, sample_csv, opentargets_response):
        init_db(db_conn)
        csv_targets = parse_seed_csv(io.StringIO(sample_csv))
        api_targets = parse_opentargets_response(opentargets_response)
        insert_targets(db_conn, csv_targets, api_targets)
        cursor = db_conn.execute(
            "SELECT gene_symbol FROM keloid_targets ORDER BY gene_symbol"
        )
        symbols = [row["gene_symbol"] for row in cursor.fetchall()]
        assert symbols == ["MTOR", "NEDD4", "PHLDA3", "TGFB1"]

    def test_total_count_with_overlap(self, db_conn, sample_csv, opentargets_response):
        init_db(db_conn)
        csv_targets = parse_seed_csv(io.StringIO(sample_csv))
        api_targets = parse_opentargets_response(opentargets_response)
        insert_targets(db_conn, csv_targets, api_targets)
        cursor = db_conn.execute("SELECT COUNT(*) as cnt FROM keloid_targets")
        assert cursor.fetchone()["cnt"] == 4


def test_seed_stores_target_role(db_conn):
    """target_role from CSV is stored in keloid_targets."""
    from scripts.seed_keloid_targets import parse_seed_csv, insert_targets
    from scripts.db import migrate_directionality_columns
    import tempfile, os, csv

    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
    writer = csv.DictWriter(tmp, fieldnames=[
        'gene_symbol','target_name','pathway','evidence_type',
        'evidence_strength','source','target_role'
    ])
    writer.writeheader()
    writer.writerow({
        'gene_symbol': 'MTOR', 'target_name': 'mTOR kinase',
        'pathway': 'PI3K/AKT/mTOR', 'evidence_type': 'functional',
        'evidence_strength': '0.7', 'source': 'manual_curation',
        'target_role': 'pro_keloid'
    })
    tmp.close()

    init_db(db_conn)
    migrate_directionality_columns(db_conn)
    with open(tmp.name) as f:
        targets = parse_seed_csv(f)
    insert_targets(db_conn, targets, [])

    row = db_conn.execute(
        "SELECT target_role FROM keloid_targets WHERE gene_symbol = 'MTOR'"
    ).fetchone()
    os.unlink(tmp.name)
    assert row is not None
    assert row['target_role'] == 'pro_keloid'
