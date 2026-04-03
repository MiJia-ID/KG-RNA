from typing import List, Dict
from Bio import SeqIO
from src.utils.logger import setup_logger

logger = setup_logger("fasta_parser")

def parse_fasta_with_metadata(fasta_file: str) -> List[Dict]:
    """
    解析包含 RNA-RBP-PMID 元数据的 FASTA 文件
    
    FASTA header 格式: >RNA_NAME|RBP_NAME|PMID1,PMID2,...
    例如: >MALAT1|HNRNPK|23555303,24284625
    
    Args:
        fasta_file: FASTA 文件路径
    
    Returns:
        候选列表，每项包含 rna, rbp, pmids, sequence
    """
    candidates = []
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # 解析 header: >RNA|RBP|PMIDs
            parts = record.id.split('|')
            
            if len(parts) < 2:
                logger.warning(f"跳过格式错误的条目: {record.id}")
                continue
            
            rna = parts[0].strip()
            rbp = parts[1].strip()
            
            # 提取 PMIDs（如果有）
            pmids = []
            if len(parts) >= 3:
                pmid_str = parts[2].strip()
                pmids = [p.strip() for p in pmid_str.split(',') if p.strip()]
            
            candidates.append({
                'rna': rna,
                'rbp': rbp,
                'pmids': pmids,
                'rna_sequence': str(record.seq)
            })
            
            logger.info(f"解析: {rna} + {rbp} ({len(pmids)} 篇文献, 序列长度: {len(record.seq)})")
        
        logger.info(f"从 FASTA 文件解析了 {len(candidates)} 个候选")
        return candidates
        
    except Exception as e:
        logger.error(f"解析 FASTA 文件失败: {e}")
        return []

def parse_fasta_simple(fasta_file: str) -> List[Dict]:
    """
    解析简单的 FASTA 文件（仅包含 RNA 名称）
    
    FASTA header 格式: >RNA_NAME
    
    Args:
        fasta_file: FASTA 文件路径
    
    Returns:
        候选列表，每项包含 rna 和 sequence
    """
    candidates = []
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            candidates.append({
                'rna': record.id.strip(),
                'rna_sequence': str(record.seq)
            })
            
            logger.info(f"解析 RNA: {record.id} (序列长度: {len(record.seq)})")
        
        logger.info(f"从 FASTA 文件解析了 {len(candidates)} 个 RNA")
        return candidates
        
    except Exception as e:
        logger.error(f"解析 FASTA 文件失败: {e}")
        return []