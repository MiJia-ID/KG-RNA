import pytest
from src.agent.retriever import detect_methods, fetch_article_by_pmid

def test_detect_methods():
    text = "We performed CLIP-seq and RIP experiments"
    methods = detect_methods(text)
    assert "CLIP" in methods
    assert "RIP" in methods

def test_fetch_article():
    article = fetch_article_by_pmid("31452104")
    assert article.get("pmid") == "31452104"
    assert "title" in article