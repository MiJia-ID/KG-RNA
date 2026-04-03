import pandas as pd
import sys
import re

def auto_find_column(columns, candidates):
    lowered = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lowered:
            return lowered[cand.lower()]
    return None

def clean_sequence(seq: str) -> str:
    if not isinstance(seq, str):
        return ""
    seq = seq.strip().replace(" ", "")
    # 只保留标准氨基酸字母
    seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())
    return seq

def export(csv_in: str, fasta_out: str,
           protein_col: str = None,
           sequence_col: str = None):
    df = pd.read_csv(csv_in, dtype=str, low_memory=False)
    cols = list(df.columns)

    # 自动匹配列
    if not protein_col:
        protein_col = auto_find_column(cols, [
            "Protein_ID (UniProt)", "protein", "Protein", "uniprot", "uniprot_id"
        ])
    if not sequence_col:
        sequence_col = auto_find_column(cols, [
            "sequence", "protein_sequence", "aa_sequence", "seq"
        ])

    if not protein_col or not sequence_col:
        print("无法匹配所需列")
        print("检测到的列:", cols)
        print(f"识别的 protein_col={protein_col}, sequence_col={sequence_col}")
        sys.exit(1)

    print(f"使用列: protein_col='{protein_col}' sequence_col='{sequence_col}'")

    n_written = 0
    with open(fasta_out, "w") as f:
        for _, r in df.iterrows():
            pid = str(r.get(protein_col, "")).strip()
            seq_raw = r.get(sequence_col, "")
            seq = clean_sequence(seq_raw)
            if not pid or not seq:
                continue
            if len(seq) < 30:  # 丢弃极短或可疑序列
                continue
            # FASTA 头可扩展附加 dataset 或 Relation
            rel = str(r.get("Relation", "")).strip()
            dataset = str(r.get("dataset", "")).strip()
            header = pid
            meta_parts = []
            if rel:
                meta_parts.append(f"REL={rel}")
            if dataset:
                meta_parts.append(f"DS={dataset}")
            if meta_parts:
                header = header + " " + " ".join(meta_parts)
            f.write(f">{header}\n{seq}\n")
            n_written += 1

    print(f"导出完成: {fasta_out} 写入序列数={n_written}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python export_kg_proteins_fasta.py <kg.csv> <out.fasta> [protein_col] [sequence_col]")
        sys.exit(1)
    csv_in = sys.argv[1]
    fasta_out = sys.argv[2]
    protein_col = sys.argv[3] if len(sys.argv) > 3 else None
    sequence_col = sys.argv[4] if len(sys.argv) > 4 else None
    export(csv_in, fasta_out, protein_col, sequence_col)
