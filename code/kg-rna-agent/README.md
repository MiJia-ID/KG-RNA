# KG-RNA Agent

基于知识图谱的 RNA-RBP 结合可行性智能评估系统

## 功能特性

- 🔍 **外部知识检索**: 从 PubMed/EuropePMC 自动检索相关文献证据
- 🤖 **LLM 推理**: 利用大语言模型进行生物学可行性推理
- 📊 **实验可信度打分**: 基于实验方法类型的加权评分系统
- 🎯 **融合评分**: 结合文献证据可信度与 LLM 推理的综合评分
- 📝 **可解释性**: 为每个候选提供详细的推理路径和证据链

## 快速开始

```bash
# 安装依赖
pip install -r requirements.txt

# 设置环境变量
export OPENAI_API_KEY="your_openai_api_key"

# 运行评估
python run_agent.py
```

## 项目结构

```
kg-rna-agent/
├── config.yaml          # 配置文件
├── requirements.txt     # 依赖
├── run_agent.py        # 主入口
├── src/
│   ├── agent/          # Agent 核心模块
│   │   ├── agent.py
│   │   ├── retriever.py
│   │   ├── llm_client.py
│   │   └── scorer.py
│   └── utils/          # 工具模块
│       └── logger.py
├── tests/              # 测试
└── outputs/            # 输出结果
```

## 输出说明

评估结果包含：
- `final_score`: 综合得分 (0-1)
- `evidence_confidence`: 实验证据可信度
- `feasibility_score`: LLM 推理可行性评分
- `llm_reasoning`: LLM 的推理结论和推荐实验