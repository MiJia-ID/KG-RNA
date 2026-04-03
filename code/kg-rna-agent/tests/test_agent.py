import pytest
from src.agent import RNAAgent

def test_agent_initialization():
    agent = RNAAgent()
    assert agent.scorer is not None

def test_evaluate_candidate_no_pmids():
    agent = RNAAgent()
    result = agent.evaluate_candidate("TEST_RNA", "TEST_RBP", [])
    assert result["rna"] == "TEST_RNA"
    assert result["rbp"] == "TEST_RBP"
    assert "final_score" in result
