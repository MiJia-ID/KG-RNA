import subprocess
import shutil
import os
from typing import Optional
from src.utils.logger import setup_logger

logger = setup_logger("diamond")

def ensure_diamond():
    exe = shutil.which("diamond")
    if not exe:
        logger.error("未找到 diamond 可执行文件，请先安装: conda install -c bioconda diamond")
        raise RuntimeError("diamond not installed")
    return exe

def build_db(fasta_in: str, db_out: str):
    ensure_diamond()
    if not os.path.exists(fasta_in):
        raise FileNotFoundError(f"KG FASTA 不存在: {fasta_in}")
    cmd = ["diamond", "makedb", "--in", fasta_in, "-d", db_out]
    logger.info(f"构建 DIAMOND 数据库: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    logger.info("数据库构建完成")

def run_blastp(query_fasta: str,
               db_path: str,
               output_file: str = "outputs/diamond_hits.tsv",
               max_target_seqs: int = 5,
               threads: int = 4,
               evalue: float = 1e-5):
    """
    运行 DIAMOND，比对失败时抛出异常，成功返回输出文件路径
    """
    # 检查数据库
    if not os.path.exists(f"{db_path}.dmnd"):
        raise FileNotFoundError(f"DIAMOND 数据库不存在: {db_path}.dmnd\n请先执行: diamond makedb --in {db_path}.fasta -d {db_path}")
    # 检查查询文件
    if not os.path.exists(query_fasta):
        raise FileNotFoundError(f"查询文件不存在: {query_fasta}")
    # 输出目录
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)

    cmd = [
        "diamond", "blastp",
        "-d", db_path,
        "-q", query_fasta,
        "-o", output_file,
        "--outfmt", "6",
        "--max-target-seqs", str(max_target_seqs),
        "--threads", str(threads),
        "--evalue", str(evalue)
    ]
    logger.info("执行: " + " ".join(cmd))
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if res.stdout:
            logger.debug(res.stdout[:500])
        if not os.path.exists(output_file):
            raise RuntimeError("比对完成但结果文件缺失: " + output_file)
        logger.info(f"DIAMOND 完成: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        logger.error(f"DIAMOND 失败，退出码={e.returncode}")
        if e.stderr:
            logger.error(e.stderr[:800])
        raise RuntimeError("DIAMOND 执行失败，请检查数据库与输入 FASTA 格式")