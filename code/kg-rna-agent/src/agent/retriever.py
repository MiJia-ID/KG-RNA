import requests
from typing import List, Dict
from src.utils.logger import setup_logger

logger = setup_logger("retriever")

EPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

METHOD_KEYWORDS = {
    "CLIP": ["CLIP", "PAR-CLIP", "ECLIP", "iCLIP", "HITS-CLIP", "eCLIP"],
    "RIP": ["RIP", "RIP-SEQ", "RIP-Seq"],
    "COIP": ["CO-IP", "COIP", "CO IMMUNOPRECIPITATION", "COIMMUNOPRECIPITATE"],
    "EMSA": ["EMSA", "GEL SHIFT", "ELECTROPHORETIC MOBILITY SHIFT"],
    "REPORTER": ["LUCIFERASE", "REPORTER ASSAY"],
    "IP": ["IMMUNOPRECIPITATION"],
    "COMPUTATIONAL": ["PREDICT", "IN SILICO", "MOTIF", "BINDING MOTIF"],
    "COPRECIP": ["CO-PRECIPITATION"]
}

def fetch_article_by_pmid(pmid: str, timeout: int = 15) -> Dict:
    """
    根据 PMID 从 EuropePMC 获取文献详情
    """
    params = {
        "query": f"EXT_ID:{pmid} OR PMID:{pmid}",
        "format": "json",
        "pageSize": 1
    }
    
    try:
        logger.info(f"Fetching PMID: {pmid}")
        r = requests.get(EPMC_SEARCH_URL, params=params, timeout=timeout)
        r.raise_for_status()
        data = r.json()
        
        if int(data.get("hitCount", 0)) == 0:
            logger.warning(f"PMID {pmid} not found in EuropePMC")
            return {}
        
        rec = data["resultList"]["result"][0]
        article = {
            "pmid": pmid,
            "title": rec.get("title", ""),
            "abstract": rec.get("abstractText", ""),
            "journal": rec.get("journalTitle", ""),
            "year": rec.get("pubYear", "")
        }
        
        logger.info(f"Successfully fetched PMID {pmid}: {article['title'][:50]}...")
        if article['abstract']:
            logger.debug(f"Abstract length: {len(article['abstract'])} chars")
        else:
            logger.warning(f"No abstract found for PMID {pmid}")
        
        return article
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to fetch PMID {pmid}: {type(e).__name__}: {e}")
        return {}

def detect_methods(text: str) -> List[str]:
    """
    从文本中检测实验方法关键词
    """
    if not text:
        return []
    
    found = set()
    t = text.upper()
    
    for method, keywords in METHOD_KEYWORDS.items():
        for kw in keywords:
            if kw.upper() in t:
                found.add(method)
                logger.debug(f"Detected method '{method}' via keyword '{kw}'")
                break
    
    return list(found)

def fetch_evidence_for_pmids(pmids: List[str]) -> List[Dict]:
    """
    批量获取多个 PMID 的证据并检测实验方法
    """
    logger.info(f"Fetching evidence for {len(pmids)} PMIDs: {pmids}")
    evidence_list = []
    
    for pmid in pmids:
        article = fetch_article_by_pmid(pmid)
        if not article:
            logger.warning(f"Skipping PMID {pmid} - no data retrieved")
            continue
        
        full_text = f"{article.get('title', '')} {article.get('abstract', '')}"
        article["detected_methods"] = detect_methods(full_text)
        
        evidence_list.append(article)
        logger.info(f"Retrieved PMID {pmid}, detected methods: {article['detected_methods'] or 'None'}")
    
    logger.info(f"Total evidence collected: {len(evidence_list)} articles")
    return evidence_list