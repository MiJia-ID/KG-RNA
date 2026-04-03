from src.agent.retriever import fetch_evidence_for_pmids

# 使用一些真实的 RNA-RBP 研究 PMID
# 这些是关于 MALAT1-HNRNPK, NEAT1-FUS 等的研究
pmids = [
    "23555303",  # MALAT1 相关
    "25599403",  # lncRNA-RBP 相关  
    "24284625"   # NEAT1 相关
]

print("=== 测试文献检索 ===\n")
evidence = fetch_evidence_for_pmids(pmids)

for ev in evidence:
    print(f"\n{'='*60}")
    print(f"PMID: {ev['pmid']}")
    print(f"Year: {ev.get('year', 'N/A')}")
    print(f"Title: {ev['title']}")
    print(f"Detected Methods: {ev['detected_methods']}")
    if ev['abstract']:
        print(f"Abstract (first 300 chars): {ev['abstract'][:300]}...")
    else:
        print("Abstract: NOT FOUND")
    print(f"{'='*60}")

print(f"\n总共检索到 {len(evidence)} 篇文献")