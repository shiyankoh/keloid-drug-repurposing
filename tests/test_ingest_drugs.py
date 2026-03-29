import json
import os
import pytest
from scripts.db import init_db
from scripts.ingest_drugs import (
    parse_dgidb_response,
    insert_drugs_and_targets,
)


@pytest.fixture
def dgidb_response():
    fixture_path = os.path.join(
        os.path.dirname(__file__), "fixtures", "dgidb_interactions.json"
    )
    with open(fixture_path) as f:
        return json.load(f)


class TestParseDgidbResponse:
    def test_filters_to_approved_only(self, dgidb_response):
        drugs, interactions = parse_dgidb_response(dgidb_response)
        drug_names = [d["drug_name"] for d in drugs]
        assert "VISTUSERTIB" not in drug_names
        assert "EVEROLIMUS" in drug_names
        assert "SIROLIMUS" in drug_names

    def test_deduplicates_drugs(self, dgidb_response):
        drugs, interactions = parse_dgidb_response(dgidb_response)
        everolimus_count = sum(1 for d in drugs if d["drug_name"] == "EVEROLIMUS")
        assert everolimus_count == 1

    def test_preserves_all_interactions(self, dgidb_response):
        drugs, interactions = parse_dgidb_response(dgidb_response)
        everolimus_genes = [
            i["gene_symbol"] for i in interactions if i["drug_name"] == "EVEROLIMUS"
        ]
        assert "MTOR" in everolimus_genes
        assert "EGFR" in everolimus_genes

    def test_extracts_action_type(self, dgidb_response):
        drugs, interactions = parse_dgidb_response(dgidb_response)
        sirolimus_int = [i for i in interactions if i["drug_name"] == "SIROLIMUS"][0]
        assert sirolimus_int["action_type"] == "inhibitor"

    def test_empty_action_type_defaults_to_unknown(self, dgidb_response):
        drugs, interactions = parse_dgidb_response(dgidb_response)
        ev_egfr = [
            i for i in interactions
            if i["drug_name"] == "EVEROLIMUS" and i["gene_symbol"] == "EGFR"
        ][0]
        assert ev_egfr["action_type"] == "unknown"

    def test_empty_response_returns_empty(self):
        empty = {"data": {"genes": {"nodes": []}}}
        drugs, interactions = parse_dgidb_response(empty)
        assert drugs == []
        assert interactions == []


class TestInsertDrugsAndTargets:
    def test_inserts_drugs(self, db_conn, dgidb_response):
        init_db(db_conn)
        drugs, interactions = parse_dgidb_response(dgidb_response)
        insert_drugs_and_targets(db_conn, drugs, interactions)
        cursor = db_conn.execute("SELECT COUNT(*) as cnt FROM drugs")
        assert cursor.fetchone()["cnt"] == 3

    def test_inserts_drug_targets(self, db_conn, dgidb_response):
        init_db(db_conn)
        drugs, interactions = parse_dgidb_response(dgidb_response)
        insert_drugs_and_targets(db_conn, drugs, interactions)
        cursor = db_conn.execute("SELECT COUNT(*) as cnt FROM drug_targets")
        assert cursor.fetchone()["cnt"] == 4

    def test_drug_target_references_correct_drug_id(self, db_conn, dgidb_response):
        init_db(db_conn)
        drugs, interactions = parse_dgidb_response(dgidb_response)
        insert_drugs_and_targets(db_conn, drugs, interactions)
        cursor = db_conn.execute("""
            SELECT d.drug_name, dt.gene_symbol
            FROM drug_targets dt JOIN drugs d ON dt.drug_id = d.id
            WHERE d.drug_name = 'EVEROLIMUS'
            ORDER BY dt.gene_symbol
        """)
        rows = cursor.fetchall()
        genes = [r["gene_symbol"] for r in rows]
        assert genes == ["EGFR", "MTOR"]
