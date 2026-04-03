import re
from src.utils.logger import setup_logger
import pandas as pd

def _extract_pmids(raw_value: str):
    if not raw_value:
        return []
    parts = re.split(r'[;,\s]+', str(raw_value))
    pmids = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        p = re.sub(r'^[Pp][Mm][Ii][Dd]:?', '', p)
        if p.isdigit() and 6 <= len(p) <= 9:
            pmids.append(p)
    return list(dict.fromkeys(pmids))

class KGSequenceMatcher:
    def __init__(self, kg_csv_path: str, protein_col: str, rna_col: str, pmid_col: str):
        self.logger = setup_logger("matcher")
        self.df = pd.read_csv(kg_csv_path, dtype=str)
        self.protein_col = protein_col
        self.rna_col = rna_col
        self.pmid_col = pmid_col
        if self.pmid_col not in self.df.columns:
            self.logger.warning(f"PMID 列不存在: {self.pmid_col}")

    def expand_candidates_from_hits(self, diamond_hits_file: str):
        candidates = []
        with open(diamond_hits_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                protein = parts[1]
                identity = float(parts[2])
                align_len = int(parts[3])
                # bitscore 通常在末尾
                try:
                    bitscore = float(parts[-1])
                except:
                    bitscore = 0.0
                sub = self.df[self.df[self.protein_col] == protein]
                if sub.empty:
                    continue
                for _, row in sub.iterrows():
                    rna_id = (row.get(self.rna_col, "") or "").strip()
                    pmid_raw = row.get(self.pmid_col, "")
                    pmids = _extract_pmids(pmid_raw)
                    candidates.append({
                        "rna": rna_id,
                        "rbp": protein,
                        "pmids": pmids,
                        "match_identity": identity,
                        "match_length": align_len,
                        "match_bitscore": bitscore
                    })
        self.logger.info(f"候选生成完成: {len(candidates)} (含PMID: {sum(1 for c in candidates if c['pmids'])})")
        return candidates