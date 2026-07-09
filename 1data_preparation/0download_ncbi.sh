#!/bin/bash
# ============================================================
# 0download_ncbi.sh — 下载 Gobiiformes 线粒体基因组
# 功能: 从 NCBI Nucleotide 下载 Gobiiformes 线粒体全基因组
#       通过 NCBI E-utilities API 获取（curl 直连 API）
# 输出: data/raw/gbff/ — 分片 GenBank 文件
#       data/raw/gobi_mitogenomes_all.gbff — 合并后总文件
# 依赖: curl, python3 (有 urllib 即可)
# 用法: tmux 中运行 bash 0download_ncbi.sh
# ============================================================

# ---- 配置 ----
PROJECT_DIR="/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec"
RAW_DIR="${PROJECT_DIR}/data/raw"
GBFF_DIR="${RAW_DIR}/gbff"

export https_proxy="socks5://172.31.150.102:20170"
export http_proxy="socks5://172.31.150.102:20170"

mkdir -p "${RAW_DIR}" "${GBFF_DIR}"

echo "=========================================="
echo "  Gobiiformes 线粒体基因组下载"
echo "  输出: ${GBFF_DIR}/"
echo "=========================================="

# ---- 步骤1: 查询 UID 列表 ----
echo ""
echo "[1/3] 查询 Gobiiformes 线粒体全基因组..."
ESRCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi\
?db=nuccore\
&term=Gobiiformes%5BOrganism%5D+AND+mitochondrion%5BTitle%5D+AND+complete%5BTitle%5D+AND+15000%3A20000%5BSequence+Length%5D\
&retmax=10000&retmode=json&usehistory=y"

curl -s --max-time 90 "${ESRCH_URL}" > /tmp/gobi_esearch.json

python3 -c "
import json
d = json.load(open('/tmp/gobi_esearch.json'))
count = d['esearchresult']['count']
ids = d['esearchresult']['idlist']
print(f'总记录数: {count}')
print(f'获取 ID 数: {len(ids)}')
with open('/tmp/gobi_ids.txt', 'w') as f:
    f.write(','.join(ids))
"

# ---- 步骤2: Python 分批下载（原生支持 SOCKS5） ----
echo ""
echo "[2/3] 分批下载 GenBank..."

python3 << 'PYEOF'
import urllib.request, json, os, sys, time

# 代理
proxy = os.environ['https_proxy']
proxy_handler = urllib.request.ProxyHandler({'https': proxy, 'http': proxy})
opener = urllib.request.build_opener(proxy_handler)
urllib.request.install_opener(opener)

with open('/tmp/gobi_ids.txt') as f:
    all_ids = f.read().strip().split(',')

total = len(all_ids)
gbff_dir = "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/raw/gbff"
batch_size = 100
n_batches = (total + batch_size - 1) // batch_size

print(f"  总数: {total} | 批次: {n_batches} (每批 {batch_size})")

for i, start in enumerate(range(0, total, batch_size)):
    batch_ids = all_ids[start:start+batch_size]
    batch_end = start + len(batch_ids) - 1
    out_file = f"{gbff_dir}/batch_{start}_{batch_end}.gbff"
    
    if os.path.exists(out_file):
        with open(out_file) as f:
            loci = f.read().count('LOCUS')
        if loci > 0:
            print(f"  [{i+1}/{n_batches}] ✅ batch {start}-{batch_end} ({loci} 条，已存在)")
            continue
    
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
           f"db=nuccore&id={','.join(batch_ids)}&rettype=gb&retmode=text")
    
    print(f"  [{i+1}/{n_batches}] ⏳ batch {start}-{batch_end}...", end=' ')
    sys.stdout.flush()
    
    try:
        resp = urllib.request.urlopen(url, timeout=180)
        content = resp.read().decode('utf-8')
        with open(out_file, 'w') as f:
            f.write(content)
        loci = content.count('LOCUS')
        print(f"✅ {loci} 条")
    except Exception as e:
        print(f"❌ 失败: {e}")
    
    time.sleep(0.5)
PYEOF

# ---- 步骤3: 合并与统计 ----
echo ""
echo "[3/3] 合并 & 统计..."
GBFF_COMBINED="${RAW_DIR}/gobi_mitogenomes_all.gbff"
cat "${GBFF_DIR}/batch_"*.gbff > "${GBFF_COMBINED}" 2>/dev/null

python3 "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/1data_preparation/_count_gbff.py" "${GBFF_COMBINED}"

echo ""
echo "=========================================="
echo "  ✅ 下载完成！"
echo "  总文件: ${GBFF_COMBINED}"
echo "=========================================="
