import pytest
from src.agent.scorer import Scorer, compute_evidence_confidence

def test_scorer_initialization():
    scorer = Scorer()
    weights = scorer.get_method_weights()
    assert "CLIP" in weights
    assert weights["CLIP"] > weights["COMPUTATIONAL"]

def test_compute_evidence_confidence():
    evidence = [
        {"detected_methods": ["CLIP"]},
        {"detected_methods": ["RIP"]}
    ]
    weights = {"CLIP": 0.9, "RIP": 0.8}
    conf = compute_evidence_confidence(evidence, weights)
    assert 0.8 < conf < 1.0