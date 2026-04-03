import os
import numpy as np
from Bio import pairwise2


def get_dssp(fasta_file, pdb_path, dssp_path):
    DSSP = "mkdssp"  # 使用 DSSP 可执行文件

    def process_dssp(dssp_file):
        with open(dssp_file, "r") as f:
            lines = f.readlines()

        # 找到DSSP 数据部分的起始行
        p = 0
        while p < len(lines) and not lines[p].startswith("  #  RESIDUE AA STRUCTURE"):
            p += 1

        # 解析DSSP 数据
        dssp_data = []
        seq = ""

        for i in range(p + 1, len(lines)):
            line = lines[i]
            if len(line) < 120:
                continue  # 确保行足够长，否则跳过

            aa = line[13]
            if aa in ("!", "*"):  # 忽略无效残基
                continue
            seq += aa

            # 直接存储 DSSP 行
            dssp_data.append(line.strip())

        return seq, dssp_data

    def match_dssp(seq, dssp_data, ref_seq):
        """如果 DSSP 序列与 PDB FASTA 序列不匹配，则进行比对调整"""
        alignments = pairwise2.align.globalxx(ref_seq, seq)
        ref_seq_aligned, seq_aligned = alignments[0].seqA, alignments[0].seqB

        new_dssp = []
        dssp_idx = 0
        for aa in seq_aligned:
            if aa == "-":
                new_dssp.append("")  # 用空行填充
            else:
                new_dssp.append(dssp_data[dssp_idx])
                dssp_idx += 1

        return [new_dssp[i] for i in range(len(ref_seq_aligned)) if ref_seq_aligned[i] != "-"]

    def extract_dssp(data_path, dssp_path, ID, ref_seq):
        """提取DSSP 并存储为原始格式"""
        dssp_file = os.path.join(dssp_path, f"{ID}.dssp")
        pdb_file = os.path.join(data_path, f"{ID}.pdb")

        try:
            # 运行 DSSP
            os.system(f"{DSSP} -i {pdb_file} -o {dssp_file}")

            # 解析 DSSP 文件
            dssp_seq, dssp_data = process_dssp(dssp_file)

            # 如果 DSSP 序列与 FASTA 不匹配，则对齐
            if dssp_seq != ref_seq:
                dssp_data = match_dssp(dssp_seq, dssp_data, ref_seq)

            # 重新写入匹配后的 DSSP 文件
            with open(dssp_file, "w") as f:
                f.write("\n".join(dssp_data) + "\n")

        except Exception as e:
            print(f"Error processing {ID}: {e}")
            return None

    # 读取 FASTA 文件
    pdbfasta = {}
    with open(fasta_file) as f:
        fasta_lines = f.readlines()

    for i in range(len(fasta_lines)):
        if fasta_lines[i].startswith(">"):
            name = fasta_lines[i].strip().split(">")[1]
            seq = fasta_lines[i + 1].strip()
            pdbfasta[name] = seq

    # 处理所有 PDB
    fault_names = []
    for name, seq in pdbfasta.items():
        if extract_dssp(pdb_path, dssp_path, name, seq) is None:
            fault_names.append(name)

    return fault_names


import argparse

def main():
    parser = argparse.ArgumentParser(description="Generate features from protein sequences.")

    parser.add_argument("--fasta_file", type=str, default='/home/mijia/egpdi_evolla/af3_data/SJT1.fasta')
    parser.add_argument('--pdb_dir', type=str, default='/home/mijia/egpdi_evolla/af3_data/af3_pdb_dir/')
    parser.add_argument("--dssp_output_path", type=str, default='/home/mijia/egpdi_evolla/af3_data/SS')

    args = parser.parse_args()

    # 调用之前定义的函数生成特征

    get_dssp(args.fasta_file, args.pdb_dir, args.dssp_output_path)


if __name__ == "__main__":
    main()
