#!/usr/bin/env python3
"""
2prepare_input.py — 准备 NeMu 管道输入：构建 BLAST 库 + 验证数据完整性

功能:
  1. 为 13 个 mt 基因构建 BLAST V5 数据库（从已提取的核苷酸 CDS FASTA）
  2. 验证每个基因的输入文件完整性
  3. 输出数据统计报告

输入: data/ncbi_db/{Gene}.fasta（CDS 核苷酸序列，来自 1extract_cds）
输出: data/ncbi_db/{Gene}_BLAST（BLAST 数据库文件）
依赖: makeblastdb, seqkit
"""
import os, subprocess, sys
from Bio import SeqIO
import pandas as pd

# ---- 配置 ----
PROJECT_DIR = "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec"
NCBI_DB_DIR = f"{PROJECT_DIR}/data/ncbi_db"
NEMU_INPUT_DIR = f"{PROJECT_DIR}/data/nemu_input"
GENES = ['A6', 'A8', 'CO1', 'CO2', 'CO3', 'Cytb',
         'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

print("=" * 50)
print("  NeMu 输入准备 — BLAST 库构建 + 验证")
print("=" * 50)

# ---- Step 1: 验证输入 FASTA ----
print("\n[1/3] 验证 CDS FASTA 文件...")
report = []
for g in GENES:
    fasta = f"{NCBI_DB_DIR}/{g}.fasta"
    if not os.path.exists(fasta):
        print(f"  ❌ {g}: 文件不存在!")
        continue
    records = list(SeqIO.parse(fasta, 'fasta'))
    lengths = [len(r.seq) for r in records]
    report.append({
        'gene': g,
        'n_seqs': len(records),
        'mean_len': f"{sum(lengths)/len(lengths):.0f}" if lengths else 0,
        'min_len': min(lengths) if lengths else 0,
        'max_len': max(lengths) if lengths else 0,
    })
    print(f"  ✅ {g}: {len(records)} 条序列")

df = pd.DataFrame(report)
print(f"\n  摘要:\n{df.to_string(index=False)}")

# ---- Step 2: 构建 BLAST DB ----
print("\n[2/3] 构建 BLAST V5 数据库...")
for g in GENES:
    fasta = f"{NCBI_DB_DIR}/{g}.fasta"
    out = f"{NCBI_DB_DIR}/{g}_BLAST"
    
    if not os.path.exists(fasta):
        continue
    
    # 检查是否已构建
    if os.path.exists(f"{out}.nsq"):
        print(f"  ✅ {g}: BLAST 库已存在")
        continue
    
    cmd = f"makeblastdb -in {fasta} -dbtype nucl -out {out} -parse_seqids"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  ✅ {g}: BLAST 库构建成功")
    else:
        print(f"  ❌ {g}: 构建失败 — {result.stderr.strip()[:100]}")

# ---- Step 3: 综合报告 ----
print("\n[3/3] 综合报告")
print(f"\n{'='*50}")
print(f"  数据准备完成!")
print(f"{'='*50}")
print(f"  BLAST 数据库: {NCBI_DB_DIR}/{{Gene}}_BLAST")
print(f"  NeMu 查询蛋白: {NEMU_INPUT_DIR}/{{Gene}}/{{Gene}}__{{Species}}.fasta")
print(f"  配置检查:")
print(f"    batch_nemu.nf params.input_dir = {NEMU_INPUT_DIR}")
print(f"    batch_nemu.nf params.db_base   = {NCBI_DB_DIR}")
print(f"\n  各基因数据量:")
for g in GENES:
    n_seq = df[df['gene']==g]['n_seqs'].values
    n_input = len(os.listdir(f"{NEMU_INPUT_DIR}/{g}")) if os.path.exists(f"{NEMU_INPUT_DIR}/{g}") else 0
    print(f"    {g}: BLAST={n_seq[0] if len(n_seq)>0 else '?'}库, NeMu={n_input}查询")
print(f"{'='*50}")
