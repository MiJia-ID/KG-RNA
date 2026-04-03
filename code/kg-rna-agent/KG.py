#!/usr/bin/env python3
# 使用: python3 /home/mijia/KG-RNA/code/process/KG.py --csv /home/mijia/KG-RNA/KG-data/protein_rna_with_fasta_clean.csv --seq-file /home/mijia/KG-RNA/KG-data/test.fasta --top 5 --method diamond
# 或   python3 scripts/seq_search.py --csv data.csv --seq "MTEYK..." --top 3

import argparse
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
from tqdm import tqdm
import requests
import sys
import os
import re
import subprocess
import json
import tempfile
import shlex

def normalize_columns(df):
    rename_map = {}
    for c in df.columns:
        lc = c.lower()
        if 'protein' in lc and ('uniprot' in lc or 'protein_id' in lc or 'protein id' in lc):
            rename_map[c] = 'Protein_ID'
        if 'rna' in lc and ('rnacentral' in lc or 'mirbase' in lc or 'rna_id' in lc or 'rna id' in lc):
            rename_map[c] = 'RNA_ID'
        if 'relation' in lc:
            rename_map[c] = 'Relation'
        if 'evidence' in lc:
            rename_map[c] = 'Evidence'
        if 'dataset' in lc:
            rename_map[c] = 'dataset'
        if 'fasta_header' in lc or 'header' in lc:
            rename_map[c] = 'fasta_header'
        if 'sequence' in lc:
            rename_map[c] = 'sequence'
    df = df.rename(columns=rename_map)
    expected = {'Protein_ID','RNA_ID','Relation','Evidence','dataset','fasta_header','sequence'}
    if not expected.issubset(set(df.columns)):
        missing = expected - set(df.columns)
        raise ValueError(f"CSV 缺少列: {missing}")
    return df

def load_csv(path):
    df = pd.read_csv(path, sep=None, engine='python', dtype=str).fillna('')
    try:
        df = normalize_columns(df)
    except ValueError:
        expected = {'Protein_ID','RNA_ID','Relation','Evidence','dataset','fasta_header','sequence'}
        if not expected.issubset(set(df.columns)):
            raise
    return df

def read_input_sequence(seq_arg, seq_file):
    if seq_arg:
        return seq_arg.strip().replace('\\n','')
    if seq_file:
        records = list(SeqIO.parse(seq_file, 'fasta'))
        if not records:
            raise ValueError("FASTA 文件无记录")
        return str(records[0].seq)
    raise ValueError("必须提供 --seq 或 --seq-file")

def pairwise_identity(s1, s2):
    alns = pairwise2.align.globalxx(s1, s2, one_alignment_only=True)
    if not alns:
        return 0.0
    a1, a2, score, start, end = alns[0]
    matches = sum(1 for x,y in zip(a1, a2) if x == y and x != '-')
    aln_len = sum(1 for x,y in zip(a1, a2) if x != '-' and y != '-')
    if aln_len == 0:
        return 0.0
    return matches / aln_len

def build_db_links(prot_id, rna_id):
    links = {}
    if prot_id:
        links['UniProt'] = f"https://rest.uniprot.org/uniprotkb/{prot_id}" if re.match(r'^[A-NR-Z0-9][0-9A-Z]{5}$', prot_id, re.I) else f"https://www.uniprot.org/uniprot/{prot_id}"
    if rna_id:
        if rna_id.startswith('URS'):
            links['RNAcentral'] = f"https://rnacentral.org/entry/{rna_id}"
        elif re.search(r'miR', rna_id, re.I) or re.match(r'hsa-|mmu-|dre-', rna_id):
            links['miRBase (search)'] = f"https://www.mirbase.org/search?q={rna_id}"
            links['RNAcentral (search)'] = f"https://rnacentral.org/search?q={rna_id}"
        else:
            links['RNAcentral (search)'] = f"https://rnacentral.org/search?q={rna_id}"
    return links

def find_top_matches_pairwise(input_seq, df, top=5):
    scores = []
    for i, row in tqdm(df.iterrows(), total=len(df), desc='比对序列 (pairwise)'):
        seq = row['sequence']
        if not seq:
            scores.append((i, 0.0))
            continue
        try:
            idt = pairwise_identity(input_seq, seq)
        except Exception:
            idt = 0.0
        scores.append((i, idt))
    scores.sort(key=lambda x: x[1], reverse=True)
    top_hits = scores[:top]
    return [(df.loc[idx], score) for idx, score in top_hits]

# --- DIAMOND 相关函数 ---
def dedupe_sequences(df, out_fasta, mapping_tsv):
    seq_to_id = {}
    seqid_rows = {}
    records = []
    for idx, seq in enumerate(df['sequence'].astype(str)):
        s = seq.strip()
        if s == '':
            continue
        if s not in seq_to_id:
            sid = f"seq{len(seq_to_id)+1}"
            seq_to_id[s] = sid
            seqid_rows[sid] = []
            records.append((sid, s))
        sid = seq_to_id[s]
        seqid_rows[sid].append(str(idx))
    with open(out_fasta, 'w') as fh:
        for sid, s in records:
            fh.write(f">{sid}\n")
            fh.write(s + "\n")
    mapping = []
    for sid, idxs in seqid_rows.items():
        mapping.append({'seqid': sid, 'indices': ';'.join(idxs)})
    pd.DataFrame(mapping).to_csv(mapping_tsv, sep='\t', index=False)
    return out_fasta, mapping_tsv

def _run_subproc(cmd_list, use_shell=False):
    if use_shell:
        return subprocess.run(cmd_list, shell=True, check=True, executable='/bin/bash')
    return subprocess.run(cmd_list, check=True)

def _diamond_available():
    try:
        subprocess.run(['diamond','version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except Exception:
        return False

def _conda_available():
    try:
        subprocess.run(['conda','--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except Exception:
        return False

def run_diamond_cmd(cmd_args, conda_env):
    """
    尝试在子进程中运行 diamond 命令：
    1) 若 diamond 可直接调用，直接执行 cmd_args 列表
    2) 否则若 conda 可用，使用 `conda run -n <env> ...`
    3) 回退：在 bash -lc 中 source conda 并 conda activate <env> 后执行命令
    """
    if _diamond_available():
        return _run_subproc(cmd_args)
    if _conda_available():
        try:
            cmd = ['conda','run','-n',conda_env] + cmd_args
            return _run_subproc(cmd)
        except subprocess.CalledProcessError:
            pass
    try:
        safe_cmd = ' '.join(shlex.quote(x) for x in cmd_args)
        try:
            base = subprocess.check_output(['conda','info','--base'], text=True).strip()
            source_cmd = f"source {shlex.quote(base)}/etc/profile.d/conda.sh && conda activate {shlex.quote(conda_env)} && {safe_cmd}"
        except Exception:
            source_cmd = f"conda activate {shlex.quote(conda_env)} && {safe_cmd}"
        return _run_subproc(source_cmd, use_shell=True)
    except Exception as e:
        raise RuntimeError(f"无法运行 diamond：{e}")

def ensure_diamond_installed(conda_env):
    if _diamond_available():
        return
    if not _conda_available():
        raise RuntimeError("diamond 不在 PATH 且 conda 不可用，请安装 diamond 或确保 conda 可用")
    try:
        subprocess.run(['conda','run','-n',conda_env,'diamond','version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return
    except Exception:
        try:
            base = subprocess.check_output(['conda','info','--base'], text=True).strip()
            cmd = f"source {shlex.quote(base)}/etc/profile.d/conda.sh && conda activate {shlex.quote(conda_env)} && diamond version"
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            return
        except Exception:
            raise RuntimeError(f"diamond 未找到，且无法在 conda 环境 {conda_env} 中执行 diamond。请安装 diamond 或确认环境名正确。")

def build_diamond_db(fasta, db_prefix, force=False, conda_env='diamond'):
    db_file = db_prefix + '.dmnd'
    if os.path.exists(db_file) and not force:
        return db_prefix
    ensure_diamond_installed(conda_env)
    cmd_args = ['diamond','makedb','--in',fasta,'-d',db_prefix]
    run_diamond_cmd(cmd_args, conda_env)
    return db_prefix

def run_diamond_search(db_prefix, query_fasta, out_tsv, top=5, threads=4, conda_env='diamond', min_id=10, min_cov=10):
    cmd_args = [
        'diamond', 'blastp',
        '--db', db_prefix,
        '--query', query_fasta,
        '--out', out_tsv,
        '--outfmt', '6',
        '--max-target-seqs', str(top),
        '--threads', str(threads),
        '--more-sensitive',
        '--id', str(min_id),
        '--query-cover', str(min_cov),
        '--subject-cover', str(min_cov),
        '--matrix', 'BLOSUM45'
    ]
    run_diamond_cmd(cmd_args, conda_env)
    return out_tsv

def parse_diamond_out(out_tsv, mapping_tsv, df, topN):
    hits = []
    map_df = pd.read_csv(mapping_tsv, sep='\t', dtype=str)
    seqid_to_indices = {r['seqid']: r['indices'].split(';') for _, r in map_df.iterrows()}
    if not os.path.exists(out_tsv):
        return hits
    # DIAMOND 默认 outfmt 6 列（12 列，按 BLAST tabular 标准顺序）
    d = pd.read_csv(out_tsv, sep='\t', header=None, dtype=str)
    # 列索引说明（BLAST tabular default）：
    # 0:qseqid 1:sseqid 2:pident 3:length ... 10:evalue 11:bitscore
    for _, row in d.iterrows():
        s = row.iloc[1]
        try:
            pident = float(row.iloc[2]) / 100.0
        except Exception:
            pident = 0.0
        bitscore = row.iloc[11] if len(row) > 11 else ''
        evalue = row.iloc[10] if len(row) > 10 else ''
        if s not in seqid_to_indices:
            continue
        for idx in seqid_to_indices[s]:
            orig_idx = int(idx)
            orig_row = df.iloc[orig_idx].to_dict()
            hit = {'orig_index': orig_idx, 'Protein_ID': orig_row.get('Protein_ID',''),
                   'RNA_ID': orig_row.get('RNA_ID',''),
                   'Relation': orig_row.get('Relation',''),
                   'Evidence': orig_row.get('Evidence',''),
                   'dataset': orig_row.get('dataset',''),
                   'fasta_header': orig_row.get('fasta_header',''),
                   'sequence': orig_row.get('sequence',''),
                   'pident': round(pident*100,3),
                   'bitscore': bitscore,
                   'evalue': evalue}
            hits.append(hit)
    df_hits = pd.DataFrame(hits)
    if df_hits.empty:
        return []
    df_hits = df_hits.sort_values(['bitscore','pident'], ascending=[False,False])
    df_hits = df_hits.drop_duplicates(subset=['orig_index'])
    return df_hits.head(topN).to_dict('records')

# --- 汇总与输出 ---
def summarize_hit(row, identity, df):
    prot = row.get('Protein_ID','')
    rna = row.get('RNA_ID','')
    info = {
        'Protein_ID': prot,
        'RNA_ID': rna,
        'identity': round(identity*100, 3) if isinstance(identity, float) else identity,
        'Relation': row.get('Relation',''),
        'Evidence': row.get('Evidence',''),
        'dataset': row.get('dataset',''),
        'fasta_header': row.get('fasta_header',''),
        'sequence': row.get('sequence',''),
        'links': build_db_links(prot, rna)
    }
    prot_rows = df[df['Protein_ID'] == prot]
    rna_list = prot_rows[['RNA_ID','Relation','Evidence','dataset']].drop_duplicates().to_dict('records')
    info['all_RNAs_for_protein'] = rna_list
    co_binders = {}
    for entry in rna_list:
        rid = entry['RNA_ID']
        others = df[(df['RNA_ID'] == rid) & (df['Protein_ID'] != prot)]['Protein_ID'].unique().tolist()
        co_binders[rid] = others
    info['co_binders'] = co_binders
    return info

def print_result(hit_infos):
    print(json.dumps(hit_infos, ensure_ascii=False, indent=2))

def build_indices(df):
    """预建索引，避免在 summarize 时全表扫描"""
    prot_to_rna = {}
    for prot, g in df.groupby('Protein_ID'):
        prot_to_rna[prot] = g[['RNA_ID','Relation','Evidence','dataset']].drop_duplicates().to_dict('records')
    rna_to_proteins = {}
    for rid, g in df.groupby('RNA_ID'):
        rna_to_proteins[rid] = g['Protein_ID'].dropna().unique().tolist()
    return {'prot_to_rna': prot_to_rna, 'rna_to_proteins': rna_to_proteins}

def summarize_hit_indexed(hit, indices, df):
    """使用预建索引快速汇总单条命中"""
    prot = hit.get('Protein_ID','')
    rna = hit.get('RNA_ID','')
    info = {
        'Protein_ID': prot,
        'RNA_ID': rna,
        'identity': hit.get('pident', ''),
        'Relation': hit.get('Relation',''),
        'Evidence': hit.get('Evidence',''),
        'dataset': hit.get('dataset',''),
        'fasta_header': hit.get('fasta_header',''),
        'sequence': hit.get('sequence',''),
        'links': build_db_links(prot, rna)
    }
    prot_rows = indices['prot_to_rna'].get(prot, [])
    info['all_RNAs_for_protein'] = prot_rows
    co_binders = {}
    for entry in prot_rows:
        rid = entry.get('RNA_ID','')
        others = [p for p in indices['rna_to_proteins'].get(rid, []) if p != prot]
        co_binders[rid] = others
    info['co_binders'] = co_binders
    return info

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--csv', required=True, help='CSV 文件路径')
    p.add_argument('--seq', help='作为查询的蛋白序列（单行字符串）')
    p.add_argument('--seq-file', help='查询序列的 fasta 文件（取第一条）')
    p.add_argument('--top', type=int, default=3, help='返回 top N 匹配')
    p.add_argument('--method', choices=['pairwise','diamond'], default='pairwise', help='比对方法')
    p.add_argument('--threads', type=int, default=4, help='diamond 线程数')
    p.add_argument('--workdir', default='.', help='临时/输出目录（diamond 模式）')
    p.add_argument('--force-build', action='store_true', help='强制重建 diamond 数据库')
    p.add_argument('--diamond-env', default='diamond', help='仅用于 diamond 步骤的 conda 环境名（默认 diamond）')
    p.add_argument('--min-id', type=float, default=10.0, help='DIAMOND 最小 identity 百分比（默认 10）')
    p.add_argument('--min-coverage', type=float, default=10.0, help='DIAMOND 最小 coverage 百分比（query/subject，默认 10）')
    args = p.parse_args()

    df = load_csv(args.csv)
    seq = read_input_sequence(args.seq, args.seq_file)

    if args.method == 'pairwise':
        hits = find_top_matches_pairwise(seq, df, top=args.top)
        results = []
        for row, score in hits:
            info = summarize_hit(row, score, df)
            results.append(info)
        print_result(results)
        return

    os.makedirs(args.workdir, exist_ok=True)
    if args.seq:
        qf = os.path.join(args.workdir, 'query.fasta')
        with open(qf, 'w') as fh:
            fh.write(">query\n")
            fh.write(args.seq.strip() + "\n")
    else:
        qf = args.seq_file

    dedup_fasta = os.path.join(args.workdir, 'deduped.fasta')
    mapping_tsv = os.path.join(args.workdir, 'seqid_map.tsv')
    dedupe_sequences(df, dedup_fasta, mapping_tsv)

    db_prefix = os.path.join(args.workdir, 'dedup_db')
    build_diamond_db(dedup_fasta, db_prefix, force=args.force_build, conda_env=args.diamond_env)

    out_tsv = os.path.join(args.workdir, 'diamond_out.tsv')
    run_diamond_search(db_prefix, qf, out_tsv, top=args.top, threads=args.threads, conda_env=args.diamond_env, min_id=args.min_id, min_cov=args.min_coverage)

    print("1) DIAMOND 完成，准备解析输出文件:", out_tsv)
    # 确认文件存在/大小
    if not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0:
        print("警告：DIAMOND 输出文件不存在或为空，退出。")
        return

    print("构建内存索引（加速后续汇总）...")
    indices = build_indices(df)

    print("2) 解析 DIAMOND 输出（parse_diamond_out）...")
    hit_records = parse_diamond_out(out_tsv, mapping_tsv, df, topN=args.top)
    print(f"   解析完成，命中数量: {len(hit_records)}")

    print("3) 汇总命中并生成结果（使用索引，加速）...")
    results = []
    for i, hr in enumerate(hit_records, 1):
        # 降低打印频率以避免大量 IO 导致变慢
        if i == 1 or i % 10 == 0:
            print(f"   处理命中 {i}/{len(hit_records)}: orig_index={hr.get('orig_index')}, Protein_ID={hr.get('Protein_ID')}")
        results.append(summarize_hit_indexed(hr, indices, df))

    out_json = os.path.join(args.workdir, 'search_results.json')
    print("4) 将结果写入文件:", out_json)
    with open(out_json, 'w', encoding='utf-8') as fh:
        json.dump(results, fh, ensure_ascii=False, indent=2)
    print("完成。已将结果写入文件。若需在终端查看前 N 行：")
    print(f"  head -n 200 {out_json}")

if __name__ == '__main__':
    main()