#!/usr/bin/env python3
"""
3rebuild_blast_with_taxid.py — 为 BLAST 库添加 taxid 映射后重建

功能: 为每个基因的 CDS FASTA 创建 taxid_map 文件，然后重建 BLAST V5 数据库。
      使 TBLASTN 的 -taxidlist 过滤功能可以正常工作。

输入: data/ncbi_db/{Gene}.fasta（CDS 核苷酸 FASTA，序列头格式: ACCESSION|SPECIES）
输出: data/ncbi_db/{Gene}_BLAST*（重建后的 BLAST 数据库）
依赖: taxonkit（需已安装 taxdump）
"""
import os, subprocess, sys

PROJECT_DIR = "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec"
NCBI_DB_DIR = f"{PROJECT_DIR}/data/ncbi_db"
GENES = ['A6', 'A8', 'CO1', 'CO2', 'CO3', 'Cytb',
         'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

for g in GENES:
    fasta = f"{NCBI_DB_DIR}/{g}.fasta"
    if not os.path.exists(fasta):
        print(f"  ⚠️ {g}: 文件不存在，跳过")
        continue
    
    # 创建 taxid_map: seq_id<tab>taxid
    map_file = f"{NCBI_DB_DIR}/{g}_taxid_map.txt"
    unmapped = set()
    
    with open(fasta) as f, open(map_file, 'w') as m:
        for line in f:
            if line.startswith('>'):
                # 格式: >PZ614122.1|Awaous_ocellaris
                seq_id = line[1:].strip().split()[0]  # PZ614122.1|Awaous_ocellaris
                species_part = seq_id.split('|')[1] if '|' in seq_id else seq_id
                species_name = species_part.replace('_', ' ')
                unmapped.add((seq_id, species_name))
    
    # 用 taxonkit 批量解析 taxid
    print(f"  {g}: 解析 {len(unmapped)} 个物种的 taxid...")
    species_names = list(set(s for _, s in unmapped))
    
    # 写物种名 → taxonkit
    proc = subprocess.run(
        ['taxonkit', 'name2taxid'],
        input='\n'.join(species_names),
        capture_output=True, text=True
    )
    
    name2taxid = {}
    for line in proc.stdout.strip().split('\n'):
        if '\t' in line:
            parts = line.split('\t')
            name2taxid[parts[0].strip()] = parts[1].strip()
    
    # 写 taxid_map
    with open(map_file, 'w') as m:
        for seq_id, species_name in unmapped:
            taxid = name2taxid.get(species_name, '0')
            m.write(f"{seq_id}\t{taxid}\n")
    
    # 统计
    n_found = sum(1 for _, s in unmapped if name2taxid.get(s, '0') != '0')
    n_miss = sum(1 for _, s in unmapped if name2taxid.get(s, '0') == '0')
    print(f"    ✓ 匹配: {n_found}, 缺失: {n_miss}")
    
    # 重建 BLAST DB
    out = f"{NCBI_DB_DIR}/{g}_BLAST"
    cmd = f"makeblastdb -in {fasta} -dbtype nucl -out {out} -parse_seqids -taxid_map {map_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"    ✓ BLAST 库重建成功")
    else:
        print(f"    ❌ 失败: {result.stderr.strip()[:200]}")

print("\n✅ 所有 BLAST 库已重建")
