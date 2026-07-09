#!/usr/bin/env python3
"""
1extract_cds.py — 从线粒体 GBFF 提取 13 个蛋白编码基因的 CDS

功能: 解析 Gobiiformes 线粒体全基因组 GenBank 文件，提取 13 个 mt 蛋白编码基因的 CDS 序列
      本脚本是 1extract_cds.ipynb 的 headless 可执行版本，功能完全一致

输入: data/raw/gobi_mitogenomes_all.gbff (675 条线粒体基因组)
输出:
  - data/ncbi_db/{Gene}.fasta  → 每个基因所有物种的 CDS 核苷酸序列（用于建 BLAST DB）
  - data/nemu_input/{Gene}/{Gene}__{Species}.fasta → 每个基因-物种的蛋白序列（NeMu 查询）
  - data/raw/cds_info.csv → 提取结果的元信息

依赖: biopython, pandas
遗传密码: 脊椎动物线粒体 (gencode=2)
"""
import os, re, sys
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# ============================================================
# Cell 1: 配置
# ============================================================
PROJECT_DIR = "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec"
GBFF_PATH = f"{PROJECT_DIR}/data/raw/gobi_mitogenomes_all.gbff"
NCBI_DB_DIR = f"{PROJECT_DIR}/data/ncbi_db"
NEMU_INPUT_DIR = f"{PROJECT_DIR}/data/nemu_input"
INFO_PATH = f"{PROJECT_DIR}/data/raw/cds_info.csv"

GENE_MAP = {
    'COX1': 'CO1', 'COI': 'CO1', 'COXI': 'CO1', 'cox1': 'CO1',
    'COX2': 'CO2', 'COII': 'CO2', 'COXII': 'CO2', 'cox2': 'CO2',
    'COX3': 'CO3', 'COIII': 'CO3', 'COXIII': 'CO3', 'cox3': 'CO3',
    'ATP6': 'A6', 'atp6': 'A6',
    'ATP8': 'A8', 'atp8': 'A8',
    'CYTB': 'Cytb', 'cob': 'Cytb',
    'ND1': 'ND1', 'nad1': 'ND1', 'NAD1': 'ND1',
    'ND2': 'ND2', 'nad2': 'ND2', 'NAD2': 'ND2',
    'ND3': 'ND3', 'nad3': 'ND3', 'NAD3': 'ND3',
    'ND4': 'ND4', 'nad4': 'ND4', 'NAD4': 'ND4',
    'ND4L': 'ND4L', 'nad4l': 'ND4L', 'NAD4L': 'ND4L',
    'ND5': 'ND5', 'nad5': 'ND5', 'NAD5': 'ND5',
    'ND6': 'ND6', 'nad6': 'ND6', 'NAD6': 'ND6',
}

TARGET_GENES = ['A6', 'A8', 'CO1', 'CO2', 'CO3', 'Cytb',
                'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

for g in TARGET_GENES:
    os.makedirs(f"{NEMU_INPUT_DIR}/{g}", exist_ok=True)
os.makedirs(NCBI_DB_DIR, exist_ok=True)

print(f"GBFF: {GBFF_PATH}")
print(f"目标基因: {TARGET_GENES}")

# ============================================================
# Cell 2: 解析 GBFF
# ============================================================
gene_data = defaultdict(list)
errors = []

records = list(SeqIO.parse(GBFF_PATH, 'genbank'))
print(f"\n总 GBFF 记录数: {len(records)}")

for rec in records:
    species = rec.annotations.get('organism', 'Unknown')
    acc = rec.id
    
    for feat in rec.features:
        if feat.type != 'CDS':
            continue
        gene_raw = feat.qualifiers.get('gene', [None])[0]
        if gene_raw is None:
            continue
        gene_std = GENE_MAP.get(gene_raw)
        if gene_std is None:
            continue
        
        translation = feat.qualifiers.get('translation', [None])[0]
        if translation is None:
            errors.append(f"{acc}/{species}/{gene_raw}: 无 translation")
            continue
        
        try:
            nucl_seq = str(feat.extract(rec.seq))
        except Exception as e:
            errors.append(f"{acc}/{species}/{gene_raw}: 提取失败 {e}")
            continue
        
        gene_data[gene_std].append((species, nucl_seq, translation, acc))

print(f"\n各基因 CDS 数量:")
for g in TARGET_GENES:
    print(f"  {g}: {len(gene_data.get(g, []))} 条")
if errors:
    print(f"\n解析错误: {len(errors)} 条")
    for e in errors[:5]:
        print(f"  ⚠️ {e}")

# ============================================================
# Cell 3: 去重——每个基因-物种保留最佳序列
# ============================================================
gene_best = {}
summary = []

for g in TARGET_GENES:
    by_species = defaultdict(list)
    for species, nucl, prot, acc in gene_data[g]:
        by_species[species].append((prot, nucl, acc))
    
    best = {}
    for species, entries in by_species.items():
        entries.sort(key=lambda x: len(x[0]), reverse=True)
        best_prot, best_nucl, best_acc = entries[0]
        best[species] = (best_nucl, best_prot, best_acc)
        summary.append({
            'gene': g, 'species': species,
            'nucl_len': len(best_nucl), 'prot_len': len(best_prot),
            'accession': best_acc, 'n_copies': len(entries),
        })
    gene_best[g] = best

df = pd.DataFrame(summary)
print(f"\n总基因-物种组合: {len(df)}")
print(f"\n各基因唯一物种数:")
for g in TARGET_GENES:
    print(f"  {g}: {len(gene_best[g])} 物种")

df.to_csv(INFO_PATH, index=False)
print(f"信息已保存: {INFO_PATH}")

# ============================================================
# Cell 4: 输出核苷酸 CDS（BLAST DB）
# ============================================================
for g in TARGET_GENES:
    records_out = []
    for species, (nucl, prot, acc) in gene_best[g].items():
        seq_id = f"{acc}|{species.replace(' ', '_')}"
        rec = SeqRecord(Seq(nucl), id=seq_id, description="")
        records_out.append(rec)
    out_path = f"{NCBI_DB_DIR}/{g}.fasta"
    SeqIO.write(records_out, out_path, 'fasta')
    print(f"{g}: {len(records_out)} 条 → {out_path}")

# ============================================================
# Cell 5: 输出蛋白序列（NeMu 查询输入）
# ============================================================
total_files = 0
for g in TARGET_GENES:
    for species, (nucl, prot, acc) in gene_best[g].items():
        species_fn = species.replace(' ', '_')
        rec = SeqRecord(Seq(prot), id=acc, description=f"{species}")
        out_path = f"{NEMU_INPUT_DIR}/{g}/{g}__{species_fn}.fasta"
        SeqIO.write(rec, out_path, 'fasta')
        total_files += 1

print(f"\n输出总文件数: {total_files}")

# ============================================================
# Cell 6: 质量报告
# ============================================================
print("\n" + "=" * 50)
print("  CDS 提取质量报告")
print("=" * 50)
print(f"GBFF 总记录: {len(records)}")
print(f"目标基因数: {len(TARGET_GENES)}")
print(f"总基因-物种组合: {len(df)}")
print(f"\n单基因物种数分布:")
print(df.groupby('gene')['species'].nunique().to_string())
print(f"\n上下游文件:")
print(f"  BLAST DB 输入: {NCBI_DB_DIR}/")
print(f"  NeMu 输入:     {NEMU_INPUT_DIR}/{{Gene}}/")
print(f"  CDS 元信息:    {INFO_PATH}")
print("=" * 50)
