#!/usr/bin/env python3
"""统计 GBFF 文件中的物种和基因分布"""
import re, sys
from collections import Counter

gbff_path = sys.argv[1] if len(sys.argv) > 1 else "/dev/stdin"

with open(gbff_path) as f:
    content = f.read()

records = [r for r in content.split('//\n') if r.strip()]
print(f"  总记录数: {len(records)}")

species = set()
genes = Counter()
for rec in records:
    m = re.search(r'ORGANISM\s+(.+?)\n', rec)
    if m:
        species.add(m.group(1).strip())
    for g in re.finditer(r'/gene="([^"]+)"', rec):
        genes[g.group(1)] += 1

print(f"  不同物种: {len(species)}")
print(f"  基因分布 (前 20):")
for g, c in genes.most_common(20):
    print(f"    {g}: {c}")
