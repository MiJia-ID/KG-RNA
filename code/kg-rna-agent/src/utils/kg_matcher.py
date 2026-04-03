import pandas as pd
from typing import List, Dict, Optional
from src.utils.logger import setup_logger

logger = setup_logger("kg_matcher")

class KGMatcher:
    """
    从知识图谱 CSV 中匹配 RNA-RBP 对应的 PMID
    """
    
    def __init__(self, csv_path: str):
        """
        初始化匹配器
        
        Args:
            csv_path: 知识图谱 CSV 文件路径
        """
        self.csv_path = csv_path
        self.kg_df = None
        self._load_kg()
    
    def _load_kg(self):
        """加载知识图谱 CSV"""
        try:
            self.kg_df = pd.read_csv(self.csv_path, dtype=str, low_memory=False)
            logger.info(f"加载知识图谱: {len(self.kg_df)} 条记录")
            
            # 打印列名以便调试
            logger.info(f"CSV 列名: {list(self.kg_df.columns)}")
            
        except Exception as e:
            logger.error(f"加载知识图谱失败: {e}")
            self.kg_df = None
    
    def find_pmids(self, rna: str, rbp: str) -> List[str]:
        """
        查找 RNA-RBP 对应的 PMID
        
        Args:
            rna: RNA 名称
            rbp: RBP 名称
        
        Returns:
            PMID 列表
        """
        if self.kg_df is None:
            logger.warning("知识图谱未加载")
            return []
        
        try:
            # 假设 CSV 有 'RNA', 'Protein', 'PMID' 列
            # 根据实际 CSV 结构调整列名
            rna_col = self._find_column(['rna', 'RNA', 'rna_name'])
            rbp_col = self._find_column(['protein', 'Protein', 'RBP', 'protein_name'])
            pmid_col = self._find_column(['pmid', 'PMID', 'pubmed_id'])
            
            if not all([rna_col, rbp_col, pmid_col]):
                logger.error(f"找不到必要的列: RNA列={rna_col}, RBP列={rbp_col}, PMID列={pmid_col}")
                return []
            
            # 查询匹配的行
            mask = (self.kg_df[rna_col].str.upper() == rna.upper()) & \
                   (self.kg_df[rbp_col].str.upper() == rbp.upper())
            
            matched_rows = self.kg_df[mask]
            
            if len(matched_rows) == 0:
                logger.warning(f"在知识图谱中未找到: {rna} + {rbp}")
                return []
            
            # 提取 PMIDs
            pmids = matched_rows[pmid_col].dropna().unique().tolist()
            
            # 如果 PMID 是用分隔符连接的字符串，需要拆分
            pmid_list = []
            for pmid in pmids:
                if ',' in str(pmid):
                    pmid_list.extend([p.strip() for p in str(pmid).split(',')])
                elif ';' in str(pmid):
                    pmid_list.extend([p.strip() for p in str(pmid).split(';')])
                else:
                    pmid_list.append(str(pmid).strip())
            
            pmid_list = [p for p in pmid_list if p and p.lower() != 'nan']
            
            logger.info(f"找到 {rna} + {rbp}: {len(pmid_list)} 篇文献")
            return pmid_list
            
        except Exception as e:
            logger.error(f"查找 PMID 失败: {e}")
            return []
    
    def _find_column(self, possible_names: List[str]) -> Optional[str]:
        """查找可能的列名"""
        if self.kg_df is None:
            return None
        
        for name in possible_names:
            for col in self.kg_df.columns:
                if col.lower() == name.lower():
                    return col
        return None
    
    def enrich_candidates(self, candidates: List[Dict]) -> List[Dict]:
        """
        为候选列表添加 PMID 信息
        
        Args:
            candidates: 候选列表（包含 rna 和 rbp）
        
        Returns:
            添加了 pmids 的候选列表
        """
        enriched = []
        
        for candidate in candidates:
            rna = candidate.get('rna')
            rbp = candidate.get('rbp')
            
            if not rna or not rbp:
                logger.warning(f"跳过缺少 RNA 或 RBP 的候选: {candidate}")
                continue
            
            # 如果已有 PMIDs，保留；否则从 KG 中查找
            pmids = candidate.get('pmids', [])
            if not pmids:
                pmids = self.find_pmids(rna, rbp)
            
            enriched.append({
                **candidate,
                'pmids': pmids
            })
        
        return enriched