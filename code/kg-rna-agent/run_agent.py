import os
import json
import yaml
import argparse
import re
import pandas as pd
from src.utils.logger import setup_logger
from src.utils.diamond_runner import run_blastp
from src.utils.kg_sequence_matcher import KGSequenceMatcher
from src.agent.agent import RNAAgent

logger = setup_logger("main")

def parse_pmids(raw: str):
    if not raw:
        return []
    parts = re.split(r'[;,\s]+', str(raw))
    out = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        p = re.sub(r'^[Pp][Mm][Ii][Dd]:?', '', p)
        if p.isdigit() and 6 <= len(p) <= 9:
            out.append(p)
    # 去重保持顺序
    return list(dict.fromkeys(out))

def load_candidates(input_file: str,
                    kg_csv: str,
                    diamond_db: str,
                    seq_mode: bool,
                    protein_col: str,
                    rna_col: str,
                    pmid_col: str):
    """
    加载候选:
    - 普通模式: JSON / CSV / FASTA(带元数据)
    - seq_mode=True: 纯蛋白序列 FASTA，需 DIAMOND + KG
    """
    if seq_mode:
        hits_file = run_blastp(
            query_fasta=input_file,
            db_path=diamond_db,
            output_file="outputs/diamond_hits.tsv",
            max_target_seqs=5
        )
        if not hits_file or not os.path.exists(hits_file):
            raise RuntimeError(f"DIAMOND 输出文件无效: {hits_file}")
        matcher = KGSequenceMatcher(
            kg_csv_path=kg_csv,
            protein_col=protein_col,
            rna_col=rna_col,
            pmid_col=pmid_col
        )
        candidates = matcher.expand_candidates_from_hits(hits_file)
        logger.info(f"扩展生成候选数: {len(candidates)}")
        return candidates

    # 非 seq-only：按 CSV 读取
    if input_file.endswith(".csv"):
        df = pd.read_csv(input_file, dtype=str)
        candidates = []
        for _, row in df.iterrows():
            pmids = parse_pmids(row.get(pmid_col, ""))
            candidates.append({
                "rna": str(row.get(rna_col, "")).strip(),
                "rbp": str(row.get(protein_col, "")).strip(),
                "pmids": pmids
            })
        logger.info(f"普通模式候选数: {len(candidates)}")
        return candidates

    raise ValueError("未支持的输入格式")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="输入文件路径")
    parser.add_argument("--seq-only", action="store_true", help="纯蛋白序列 FASTA 模式")
    parser.add_argument("--kg-csv", required=True, type=str, help="知识图谱 CSV")
    parser.add_argument("--diamond-db", required=True, type=str, help="DIAMOND 数据库前缀 (不含 .dmnd 可自动补)")
    parser.add_argument("--kg-protein-col", required=True, type=str, help="KG 蛋白列名")
    parser.add_argument("--kg-rna-col", required=True, type=str, help="KG RNA 列名")
    parser.add_argument("--kg-pmid-col", required=True, type=str, help="KG PMID/Evidence 列名")
    parser.add_argument("-o", "--output", default="outputs/evaluation_results.json", type=str)
    args = parser.parse_args()

    logger.info("KG-RNA Agent started")

    candidates = load_candidates(
        input_file=args.input,
        kg_csv=args.kg_csv,
        diamond_db=args.diamond_db,
        seq_mode=args.seq_only,
        protein_col=args.kg_protein_col,
        rna_col=args.kg_rna_col,
        pmid_col=args.kg_pmid_col
    )

    agent = RNAAgent(csv_path=args.kg_csv)
    results = agent.evaluate_batch(candidates)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    logger.info(f"结果已保存: {args.output}")

    # 简单打印
    for i, r in enumerate(results, 1):
        print(f"{i}. {r['rna']} + {r['rbp']} -> {r.get('final_score', 0):.3f}")

if __name__ == "__main__":
    main()