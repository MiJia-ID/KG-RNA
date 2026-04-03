import os
import re
import json
from openai import OpenAI
from src.utils.logger import setup_logger

logger = setup_logger("llm_client")

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
OPENAI_MODEL = os.getenv("OPENAI_MODEL", "gpt-3.5-turbo")
OPENAI_BASE_URL = os.getenv("OPENAI_BASE_URL", "https://api.openai-hub.com/v1")
USE_RULE_BASED = os.getenv("USE_RULE_BASED", "false").lower() == "true"

client = None
if OPENAI_API_KEY and not USE_RULE_BASED:
    try:
        client = OpenAI(
            api_key=OPENAI_API_KEY,
            base_url=OPENAI_BASE_URL
        )
        logger.info(f"OpenAI client initialized with base_url={OPENAI_BASE_URL}")
    except Exception as e:
        logger.warning(f"Failed to init OpenAI client: {e}, falling back to rule-based scoring")

def call_llm_for_candidate(evidence_text: str, match_info: dict = None, max_tokens: int = 512, temperature: float = 0.0) -> dict:
    """
    调用 LLM 或规则评分
    """
    # 优先规则评分
    if USE_RULE_BASED or not client:
        logger.info("Using rule-based scoring (LLM bypassed)")
        return generate_rule_based_score(evidence_text, match_info)
    
    try:
        logger.info(f"Calling OpenAI API with model: {OPENAI_MODEL}")
        prompt = f"""你是 RNA-蛋白质互作专家。根据以下证据评估可行性（0-1分）：

{evidence_text}

返回 JSON 格式：
{{
  "feasibility_score": 0.0-1.0,
  "evidence_confidence": 0.0-1.0,
  "summary": "简短结论",
  "key_evidence": ["关键证据1", "关键证据2"],
  "recommended_next_experiments": ["建议实验1", "建议实验2"]
}}"""

        response = client.chat.completions.create(
            model=OPENAI_MODEL,
            messages=[{"role": "user", "content": prompt}],
            max_tokens=max_tokens,
            temperature=temperature
        )
        
        # 检查响应类型
        if not hasattr(response, 'choices'):
            logger.error(f"Unexpected response type: {type(response)}")
            logger.error(f"Response content: {response}")
            raise ValueError(f"Invalid API response format: {type(response)}")
        
        # 提取内容
        content = response.choices[0].message.content.strip()
        logger.info(f"LLM raw response: {content[:200]}...")
        
        # 解析 JSON
        json_start = content.find("{")
        json_end = content.rfind("}") + 1
        
        if json_start >= 0 and json_end > json_start:
            json_str = content[json_start:json_end]
            result = json.loads(json_str)
            
            # 验证必需字段
            required_keys = ["feasibility_score", "evidence_confidence", "summary"]
            for key in required_keys:
                if key not in result:
                    logger.warning(f"Missing key '{key}' in LLM response, using default")
                    result[key] = 0.0 if "score" in key else "信息不足"
            
            # 确保列表字段存在
            result.setdefault("key_evidence", [])
            result.setdefault("recommended_next_experiments", [])
            
            return result
        else:
            raise ValueError("LLM 未返回有效 JSON")
            
    except json.JSONDecodeError as e:
        logger.error(f"JSON parsing failed: {e}")
        logger.error(f"Content was: {content}")
        return generate_rule_based_score(evidence_text, match_info)
        
    except Exception as e:
        logger.error(f"LLM call failed ({type(e).__name__}): {str(e)}")
        # 降级到规则评分
        return generate_rule_based_score(evidence_text, match_info)

def generate_rule_based_score(evidence_text: str, match_info: dict = None) -> dict:
    """
    基于序列比对质量的规则评分
    """
    score = 0.3
    confidence = 0.2
    summary_parts = []
    
    if match_info:
        identity = match_info.get("match_identity", 0)
        length = match_info.get("match_length", 0)
        bitscore = match_info.get("match_bitscore", 0)
        
        # 评分规则
        if identity >= 90:
            score = 0.9
            confidence = 0.85
            summary_parts.append(f"高同源性 ({identity:.1f}%)")
        elif identity >= 70:
            score = 0.75
            confidence = 0.7
            summary_parts.append(f"较高同源性 ({identity:.1f}%)")
        elif identity >= 50:
            score = 0.6
            confidence = 0.5
            summary_parts.append(f"中等同源性 ({identity:.1f}%)")
        else:
            score = 0.4
            confidence = 0.3
            summary_parts.append(f"低同源性 ({identity:.1f}%)")
        
        # 比对长度调整
        if length >= 100:
            score += 0.05
            summary_parts.append(f"长比对 ({length} aa)")
        elif length < 50:
            score -= 0.1
            summary_parts.append(f"短比对 ({length} aa)")
        
        # Bitscore 调整
        if bitscore >= 200:
            score += 0.05
            summary_parts.append(f"高 Bitscore ({bitscore:.1f})")
    else:
        summary_parts.append("无比对信息")
    
    # 检查文献证据
    has_pmid = bool(re.search(r'\d{7,8}', evidence_text or ""))
    if has_pmid:
        score += 0.1
        confidence += 0.15
        summary_parts.append("含文献支持")
    else:
        summary_parts.append("缺少文献验证")
    
    # 限制范围
    score = min(1.0, max(0.0, score))
    confidence = min(1.0, max(0.0, confidence))
    
    summary = "规则评分：" + "，".join(summary_parts) if summary_parts else "基于序列同源性预测"
    
    return {
        "feasibility_score": round(score, 3),
        "evidence_confidence": round(confidence, 3),
        "summary": summary,
        "key_evidence": [
            f"序列同源性: {match_info.get('match_identity', 0):.1f}%" if match_info else "序列比对",
            f"比对长度: {match_info.get('match_length', 0)} aa" if match_info and match_info.get('match_length') else "",
            f"Bitscore: {match_info.get('match_bitscore', 0):.1f}" if match_info and match_info.get('match_bitscore') else ""
        ],
        "recommended_next_experiments": [
            "CLIP-seq (交联免疫沉淀)",
            "RIP-qPCR (RNA 免疫沉淀)",
            "EMSA (凝胶迁移)",
            "双荧光素酶报告基因"
        ]
    }