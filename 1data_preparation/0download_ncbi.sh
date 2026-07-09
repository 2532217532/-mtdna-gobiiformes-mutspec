#!/bin/bash
# ============================================================
# 0download_ncbi.sh — 下载 Gobiiformes 线粒体基因组
# 功能: 从 NCBI Nucleotide 下载 Gobiiformes 线粒体全基因组
#       通过 NCBI E-utilities API 直接获取（解决 esearch/efetch 代理问题）
# 输出: data/raw/gbff/ — GenBank 格式文件，按批次分片存储
# 依赖: curl, python3
# 用法: tmux新会话中运行 bash 0download_ncbi.sh
# ============================================================
set -euo pipefail

# ---- 配置 ----
PROJECT_DIR="/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec"
RAW_DIR="${PROJECT_DIR}/data/raw"
GBFF_DIR="${RAW_DIR}/gbff"
NCBI_DB_DIR="${PROJECT_DIR}/data/ncbi_db"
NEMU_INPUT_DIR="${PROJECT_DIR}/data/nemu_input"
GENCODE=2

# 代理
export https_proxy="socks5://172.31.150.102:20170"
export http_proxy="socks5://172.31.150.102:20170"

mkdir -p "${RAW_DIR}" "${GBFF_DIR}" "${NCBI_DB_DIR}" "${NEMU_INPUT_DIR}"

echo "=========================================="
echo "  Gobiiformes 线粒体基因组下载"
echo "  输出目录: ${GBFF_DIR}"
echo "=========================================="

# ---- 步骤1: 查询 UIDs ----
echo ""
echo "[1/4] 查询 Gobiiformes 线粒体全基因组 UID 列表..."
QUERY_ENC="Gobiiformes[Organism] AND mitochondrion[Title] AND complete[Title] AND 15000:20000[Sequence Length]"
ESRCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=${QUERY_ENC}&retmax=10000&retmode=json&usehistory=y"

curl -s --max-time 60 "${ESRCH_URL}" > /tmp/gobi_esearch.json
COUNT=$(python3 -c "import json; d=json.load(open('/tmp/gobi_esearch.json')); print(d['esearchresult']['count'])")
IDS=$(python3 -c "import json; d=json.load(open('/tmp/gobi_esearch.json')); print(','.join(d['esearchresult']['idlist']))")
WEBENV=$(python3 -c "import json; d=json.load(open('/tmp/gobi_esearch.json')); print(d['esearchresult'].get('webenv',''))")
QUERY_KEY=$(python3 -c "import json; d=json.load(open('/tmp/gobi_esearch.json')); print(d['esearchresult'].get('querykey','1'))")

echo "  总记录数: ${COUNT}"

if [ -z "$IDS" ]; then
    echo "  未能获取 ID 列表，尝试用 WebEnv 分批下载..."
    HAS_WEBENV=true
else
    echo "  UID 列表长度: $(echo $IDS | tr ',' '\n' | wc -l)"
    HAS_WEBENV=false
fi

# ---- 步骤2: 下载 GenBank 格式 ----
echo ""
echo "[2/4] 下载 GenBank 文件..."

# 方法A: 用 ID 列表直下 (快)
if [ "$HAS_WEBENV" = false ]; then
    # 分批次，每批 100 条
    IFS=',' read -ra ID_ARRAY <<< "$IDS"
    TOTAL=${#ID_ARRAY[@]}
    BATCH=100
    
    for ((start=0; start<TOTAL; start+=BATCH)); do
        end=$((start + BATCH))
        [ $end -gt $TOTAL ] && end=$TOTAL
        BATCH_IDS=$(IFS=,; echo "${ID_ARRAY[*]:start:end-start}")
        OUT_FILE="${GBFF_DIR}/batch_${start}_$((end-1)).gbff"
        
        if [ -f "$OUT_FILE" ] && [ $(grep -c '^LOCUS' "$OUT_FILE" 2>/dev/null || echo 0) -gt 0 ]; then
            echo "  已存在: batch ${start}-$((end-1))，跳过"
            continue
        fi
        
        EFTCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${BATCH_IDS}&rettype=gb&retmode=text"
        echo "  下载: ${start}-$((end-1))/${TOTAL}"
        curl -s --max-time 120 "${EFTCH_URL}" > "${OUT_FILE}"
        
        # 检查是否成功
        LINES=$(wc -l < "${OUT_FILE}")
        LOCI=$(grep -c '^LOCUS' "${OUT_FILE}" 2>/dev/null || echo 0)
        echo "    → ${LOCI} 条记录 (${LINES} 行)"
        
        # 休息 1 秒避免请求过频
        sleep 1
    done
else
    # 方法B: 用 WebEnv + query_key
    for ((start=0; start<COUNT; start+=100)); do
        end=$((start + 99))
        [ $end -ge $COUNT ] && end=$((COUNT - 1))
        OUT_FILE="${GBFF_DIR}/batch_${start}_${end}.gbff"
        
        if [ -f "$OUT_FILE" ] && [ $(grep -c '^LOCUS' "$OUT_FILE" 2>/dev/null || echo 0) -gt 0 ]; then
            echo "  已存在: batch ${start}-${end}，跳过"
            continue
        fi
        
        EFTCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&WebEnv=${WEBENV}&query_key=${QUERY_KEY}&retstart=${start}&retmax=100&rettype=gb&retmode=text"
        echo "  下载: ${start}-${end}/${COUNT}"
        curl -s --max-time 120 "${EFTCH_URL}" > "${OUT_FILE}"
        
        LINES=$(wc -l < "${OUT_FILE}")
        LOCI=$(grep -c '^LOCUS' "${OUT_FILE}" 2>/dev/null || echo 0)
        echo "    → ${LOCI} 条记录 (${LINES} 行)"
        sleep 1
    done
fi

# ---- 步骤3: 合并与统计 ----
echo ""
echo "[3/4] 合并与统计..."
cat "${GBFF_DIR}/batch_"*.gbff > "${RAW_DIR}/gobi_mitogenomes_all.gbff" 2>/dev/null || true
TOTAL_RECORDS=$(grep -c '^LOCUS' "${RAW_DIR}/gobi_mitogenomes_all.gbff" 2>/dev/null || echo 0)
echo "  总 GenBank 记录: ${TOTAL_RECORDS}"

# 提取物种名概览
if [ $TOTAL_RECORDS -gt 0 ]; then
    echo ""
    echo "  物种列表 (前30):"
    grep '^ORGANISM' "${RAW_DIR}/gobi_mitogenomes_all.gbff" | sort -u | head -30
    echo "  不同物种数: $(grep '^ORGANISM' ${RAW_DIR}/gobi_mitogenomes_all.gbff | sort -u | wc -l)"
fi

# ---- 步骤4: 提取 CDS 基因信息 ----
echo ""
echo "[4/4] CDS 基因分布概览..."
if [ $TOTAL_RECORDS -gt 0 ]; then
    python3 << 'PYEOF'
import re, sys
from collections import Counter

with open("/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/raw/gobi_mitogenomes_all.gbff") as f:
    content = f.read()

# 按 // 分隔记录
records = content.split('//')
print(f"  记录数: {len(records)}")

genes = Counter()
for rec in records:
    for m in re.finditer(r'/gene="([^"]+)"', rec):
        genes[m.group(1)] += 1

print("  基因分布:")
for g, c in genes.most_common(30):
    print(f"    {g}: {c}")
PYEOF
fi

echo ""
echo "=========================================="
echo "  ✅ 下载完成！"
echo "  GenBank: ${RAW_DIR}/gobi_mitogenomes_all.gbff"
echo "  分片:    ${GBFF_DIR}/"
echo "=========================================="
