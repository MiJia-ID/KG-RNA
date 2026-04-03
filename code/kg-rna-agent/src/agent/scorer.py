import os
import pandas as pd
from typing import List, Dict
from src.utils.logger import setup_logger

logger = setup_logger("scorer")

KEYWORD_TO_METHOD = {
    "CLIP": ["CLIP", "PAR-CLIP", "ECLIP", "ICLIP", "HITS-CLIP", "eCLIP"],
    "RIP": ["RIP", "RIP-SEQ"],
    "COIP": ["CO-IP", "COIP", "CO IMMUNOPRECIPITATION", "COIMMUNOPRECIPITATE"],
    "EMSA": ["EMSA", "GEL SHIFT", "ELECTROPHORETIC MOBILITY SHIFT"],
    "REPORTER": ["LUCIFERASE", "REPORTER ASSAY"],
    "IP": ["IMMUNOPRECIPITATION"],
    "COMPUTATIONAL": ["PREDICT", "IN SILICO", "MOTIF", "BINDING MOTIF"],
    "COPRECIP": ["CO-PRECIPITATION"]
}

BASE_WEIGHTS = {
    "CLIP": 0.92,
    "RIP": 0.82,
    "COIP": 0.75,
    "EMSA": 0.72,
    "REPORTER": 0.65,
    "IP": 0.6,
    "COMPUTATIONAL": 0.4,
    "COPRECIP": 0.7
}

def infer_methods_from_text(text: str) -> List[str]:
    """从文本中推断实验方法"""
    t = (text or "").upper()
    found = set()
    for method, kwlist in KEYWORD_TO_METHOD.items():
        for kw in kwlist:
            if kw in t:
                found.add(method)
                break
    return list(found)

def build_method_weights_from_csv(csv_path: str) -> Dict[str, float]:
    """
    从 CSV 文件的 Evidence 列统计方法并调整权重
    """
    if not os.path.exists(csv_path):
        logger.warning(f"CSV not found: {csv_path}, using base weights")
        return BASE_WEIGHTS.copy()

    try:
        df = pd.read_csv(csv_path, dtype=str, low_memory=False)
    except Exception as e:
        logger.error(f"Failed to read CSV: {e}")
        return BASE_WEIGHTS.copy()

    ev_col = None
    for c in df.columns:
        if c.lower() == "evidence":
            ev_col = c
            break

    if ev_col is None:
        logger.warning("Evidence column not found")
        return BASE_WEIGHTS.copy()

    counts = {}
    for text in df[ev_col].fillna(""):
        methods = infer_methods_from_text(text)
        for m in methods:
            counts[m] = counts.get(m, 0) + 1

    if not counts:
        return BASE_WEIGHTS.copy()

    total = sum(counts.values())
    weights = {}
    for m, cnt in counts.items():
        base = BASE_WEIGHTS.get(m, 0.45)
        freq = cnt / total
        adj = base + min(0.08, freq * 0.2)
        weights[m] = round(min(0.99, adj), 3)

    for m in BASE_WEIGHTS:
        if m not in weights:
            weights[m] = BASE_WEIGHTS[m]

    logger.info(f"Method weights built from CSV")
    return weights

def compute_evidence_confidence(evidence_list: List[Dict], method_weights: Dict[str, float]) -> float:
    """
    计算证据可信度
    """
    if not evidence_list:
        return 0.0

    scores = []
    for ev in evidence_list:
        methods = ev.get("detected_methods", []) or []
        if not methods:
            scores.append(0.3)
            continue
        w = max(method_weights.get(m, 0.45) for m in methods)
        scores.append(w)

    base = sum(scores) / len(scores)
    count_bonus = min(0.1, 0.02 * len(evidence_list))
    
    return round(min(1.0, base + count_bonus), 3)

def aggregate_score(evidence_confidence: float, feasibility_score: float, alpha: float = 0.55) -> float:
    """
    融合证据可信度与可行性评分
    """
    return round(alpha * evidence_confidence + (1 - alpha) * feasibility_score, 3)

class Scorer:
    """
    打分器：融合实验证据可信度与 LLM 推理可行性
    """
    def __init__(self, csv_path: str = None):
        self.csv_path = csv_path or "/home/mijia/KG-RNA/KG-data/protein_rna_with_fasta_clean.csv"
        self.method_weights = build_method_weights_from_csv(self.csv_path)

    def get_method_weights(self) -> Dict[str, float]:
        return self.method_weights

    def score_candidate(self, evidence_list: List[Dict], feasibility_score: float, alpha: float = 0.55) -> Dict:
        ev_conf = compute_evidence_confidence(evidence_list, self.method_weights)
        final = aggregate_score(ev_conf, feasibility_score, alpha=alpha)
        
        return {
            "evidence_confidence": ev_conf,
            "feasibility_score": round(float(feasibility_score), 3),
            "final_score": final,
            "method_weights_used": self.method_weights
        }