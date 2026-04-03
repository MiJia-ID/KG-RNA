import sys
from Bio import SeqIO

def convert_fasta(input_file, output_file, mapping_file=None):
    """
    转换 FASTA 文件格式
    
    Args:
        input_file: 输入 FASTA 文件
        output_file: 输出 FASTA 文件 (标准格式)
        mapping_file: RNA-RBP-PMID 映射文件 (可选，CSV 或 TSV)
    """
    # 如果提供了映射文件，加载映射
    mapping = {}
    if mapping_file:
        import pandas as pd
        df = pd.read_csv(mapping_file, sep=None, engine='python')
        print(f"映射文件列: {list(df.columns)}")
        
        # 假设有 id, rna, rbp, pmids 列
        for _, row in df.iterrows():
            seq_id = str(row.get('id', row.get('sequence_id', '')))
            mapping[seq_id] = {
                'rna': row.get('rna', row.get('RNA', '')),
                'rbp': row.get('rbp', row.get('RBP', row.get('protein', ''))),
                'pmids': str(row.get('pmids', row.get('PMID', ''))).replace(' ', ',')
            }
    
    # 转换
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(input_file, "fasta"):
            seq_id = record.id
            
            if seq_id in mapping:
                info = mapping[seq_id]
                new_id = f"{info['rna']}|{info['rbp']}|{info['pmids']}"
            else:
                # 如果没有映射，尝试从原 ID 提取
                new_id = f"RNA_{seq_id}|UNKNOWN|"
            
            out.write(f">{new_id}\n")
            out.write(f"{record.seq}\n")
            print(f"转换: {seq_id} -> {new_id}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: python convert_fasta.py <input.fasta> <output.fasta> [mapping.csv]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    mapping_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    convert_fasta(input_file, output_file, mapping_file)
    print(f"\n转换完成: {output_file}")
