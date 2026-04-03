import json
from typing import List, Dict
from .retriever import fetch_evidence_for_pmids
from .llm_client import call_llm_for_candidate
from .scorer import Scorer
from src.utils.logger import setup_logger

logger = setup_logger("agent")

def build_evidence_text(rna: str, rbp: str, evidence: List[Dict]) -> str:
    """
    构建发送给 LLM 的证据文本
    """
    parts = [f"=== Candidate: RNA={rna}, RBP={rbp} ===\n"]
    
    for i, ev in enumerate(evidence, 1):
        pmid = ev.get('pmid', 'unknown')
        title = ev.get('title', 'No title')
        abstract = ev.get('abstract', 'No abstract')[:800]
        methods = ", ".join(ev.get("detected_methods", [])) or "Unknown"
        year = ev.get('year', 'N/A')
        journal = ev.get('journal', 'N/A')
        
        parts.append(f"[Evidence {i}]")
        parts.append(f"PMID: {pmid}")
        parts.append(f"Year: {year}")
        parts.append(f"Journal: {journal}")
        parts.append(f"Title: {title}")
        parts.append(f"Abstract: {abstract}")
        parts.append(f"Detected methods: {methods}")
        parts.append("---\n")
    
    return "\n".join(parts)

class RNAAgent:
    """
    RNA-RBP 结合可行性评估的智能 Agent
    """
    
    def __init__(self, csv_path: str = None):
        self.scorer = Scorer(csv_path=csv_path)
        self.logger = setup_logger("agent")          # 修复
        self.retriever = fetch_evidence_for_pmids    # 简单引用函数
        self.logger.info("RNAAgent initialized")

    def evaluate_candidate(self, rna: str, rbp: str, pmids: list, match_info: dict = None):
        """
        对单个候选进行完整评估
        """
        self.logger.info(f"评估: {rna} + {rbp} PMIDs={pmids}")
        evidence_list = self.retriever(pmids) if pmids else []
        if pmids and not evidence_list:
            self.logger.warning(f"PMIDs 未检索到文献: {pmids}")

        evidence_text = build_evidence_text(rna, rbp, evidence_list)
        # 附加比对信息（如有）
        if match_info:
            evidence_text += f"\n\n[Sequence Alignment]\nIdentity={match_info.get('match_identity',0):.1f}%\nLength={match_info.get('match_length',0)}\nBitscore={match_info.get('match_bitscore',0):.1f}\n"

        llm_resp = call_llm_for_candidate(evidence_text, match_info=match_info)

        feasibility = 0.0
        if isinstance(llm_resp, dict):
            try:
                feasibility = float(llm_resp.get("feasibility_score", 0.0))
            except:
                feasibility = 0.0

        scored = self.scorer.score_candidate(evidence_list, feasibility)

        result = {
            "rna": rna,
            "rbp": rbp,
            "pmids": pmids,
            "final_score": scored["final_score"],
            "evidence_confidence": scored["evidence_confidence"],
            "feasibility_score": scored["feasibility_score"],
            "llm_reasoning": {
                "summary": llm_resp.get("summary", "LLM 推理失败"),
                "key_evidence": llm_resp.get("key_evidence", []),
                "recommended_experiments": llm_resp.get("recommended_next_experiments", []),
                "raw_response": llm_resp if "error" in llm_resp else None
            },
            "evidence_details": evidence_list,
            "method_weights_used": scored["method_weights_used"]
        }
        self.logger.info(f"完成: final_score={result['final_score']:.3f}")
        return result

    def evaluate_batch(self, candidates: List[Dict]) -> List[Dict]:
        """
        批量评估候选
        """
        logger.info(f"Batch evaluation: {len(candidates)} candidates")
        results = []
        
        for i, c in enumerate(candidates, 1):
            logger.info(f"Processing {i}/{len(candidates)}")
            res = self.evaluate_candidate(c["rna"], c["rbp"], c.get("pmids", []))
            results.append(res)
        
        results.sort(key=lambda x: x["final_score"], reverse=True)
        logger.info("Batch evaluation complete")
        return results

    def explain_top_candidate(self, result: Dict) -> str:
        """
        生成可解释性报告
        """
        explanation = f"""
{'=' * 70}
RNA-RBP 结合可行性评估报告
{'=' * 70}

【候选信息】
  RNA: {result['rna']}
  RBP: {result['rbp']}
  参考文献数: {len(result['pmids'])}

【综合评分】
  最终得分: {result['final_score']:.3f}
  ├─ 实验证据可信度: {result['evidence_confidence']:.3f}
  └─ LLM 推理可行性: {result['feasibility_score']:.3f}

【LLM 推理结论】
{result['llm_reasoning']['summary']}

【关键支持证据】"""
        
        key_ev = result['llm_reasoning']['key_evidence']
        if key_ev:
            for ev in key_ev:
                explanation += f"\n  • {ev}"
        else:
            explanation += "\n  无关键证据"
        
        explanation += "\n\n【推荐后续实验】"
        rec_exp = result['llm_reasoning']['recommended_experiments']
        if rec_exp:
            for exp in rec_exp:
                explanation += f"\n  • {exp}"
        else:
            explanation += "\n  无推荐实验"
        
        explanation += f"\n\n【检测到的实验方法】"
        all_methods = set()
        for ev in result['evidence_details']:
            all_methods.update(ev.get('detected_methods', []))
        if all_methods:
            explanation += f"\n  {', '.join(sorted(all_methods))}"
        else:
            explanation += "\n  未检测到明确方法"
        
        explanation += f"\n\n{'=' * 70}\n"
        return explanation