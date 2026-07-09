#!/bin/bash
# ============================================================
# 0download_ncbi.sh — 下载 Gobiiformes 线粒体基因组
# 功能: 从 NCBI Nucleotide 下载 Gobiiformes 线粒体全基因组
#       每个物种一个 GenBank 文件，含完整 CDS 注释
# 输出: data/raw/gbff/ — GenBank 格式原始文件
# 依赖: esearch/efetch (entrez-direct), datasets
# 用法: bash 0download_ncbi.sh
# ============================================================
set -euo pipefail

# ---- 配置 ----
PROJECT_DIR="/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec"
RAW_DIR="${PROJECT_DIR}/data/raw"
NCBI_DB_DIR="${PROJECT_DIR}/data/ncbi_db"
NEMU_INPUT_DIR="${PROJECT_DIR}/data/nemu_input"
GENCODE=2  # 脊椎动物线粒体遗传密码

# 代理设置
export https_proxy="socks5://172.31.150.102:20170"
export http_proxy="socks5://172.31.150.102:20170"

# Conda 环境
CONDA_BASE="/home/zengjl/miniforge3"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate gobi_mutspec
export JAVA_HOME="${CONDA_BASE}/envs/gobi_mutspec/jdk-17.0.14+7"
export PATH="${JAVA_HOME}/bin:${PATH}"

mkdir -p "${RAW_DIR}/gbff"

echo "=========================================="
echo "  步骤1: 查询 Gobiiformes 线粒体基因组数量"
echo "=========================================="

# 查询: Gobiiformes 的线粒体全基因组 (15-20kb)
QUERY='Gobiiformes[Organism] AND mitochondrion[Title] AND complete[Title] AND 15000:20000[Sequence Length]'
esearch -db nucleotide -query "$QUERY" > /tmp/gobi_esearch_result.txt 2>/dev/null || true

COUNT=$(grep -oP '<Count>\K[^<]+' /tmp/gobi_esearch_result.txt 2>/dev/null || echo "0")
echo "线粒体全基因组数量: ${COUNT}"

# 如果有 WebEnv，用 efetch 下 GenBank 格式
WEBENV=$(grep -oP '<WebEnv>\K[^<]+' /tmp/gobi_esearch_result.txt 2>/dev/null || echo "")
QUERY_KEY=$(grep -oP '<QueryKey>\K[^<]+' /tmp/gobi_esearch_result.txt 2>/dev/null || echo "1")

if [ -n "$WEBENV" ] && [ "$COUNT" -gt 0 ]; then
    echo ""
    echo "=========================================="
    echo "  步骤2: 下载线粒体基因组 (GenBank 格式)"
    echo "  总数: ${COUNT}"
    echo "  输出: ${RAW_DIR}/gbff/"
    echo "=========================================="
    
    # 批量下载，每次 50 条
    BATCH=50
    for ((start=0; start<COUNT; start+=BATCH)); do
        end=$((start + BATCH - 1))
        [ $end -ge $((COUNT - 1)) ] && end=$((COUNT - 1))
        OUT_FILE="${RAW_DIR}/gbff/gobi_mitogenomes_${start}_${end}.gbff"
        
        echo "  下载: ${start} - ${end} → ${OUT_FILE}"
        efetch -db nucleotide -format gb -WebEnv "${WEBENV}" -query_key "${QUERY_KEY}" \
            -start $((start + 1)) -stop $((end + 1)) > "${OUT_FILE}" 2>/dev/null || \
        efetch -db nucleotide -format gb -WebEnv "${WEBENV}" -query_key "${QUERY_KEY}" \
            -start $((start + 1)) -stop $((end + 1)) > "${OUT_FILE}" 2>/dev/null || \
        echo "  ⚠️  下载失败: batch ${start}-${end}"
        
        # 显示进度
        RECORDS=$(grep -c '^LOCUS' "${OUT_FILE}" 2>/dev/null || echo "0")
        echo "     → 获取 ${RECORDS} 条记录"
    done
    
    # 合并所有 gbff 文件
    echo ""
    echo "  合并 GBFF 文件..."
    cat "${RAW_DIR}/gbff/gobi_mitogenomes_"*.gbff > "${RAW_DIR}/gobi_mitogenomes_all.gbff" 2>/dev/null || true
    TOTAL=$(grep -c '^LOCUS' "${RAW_DIR}/gobi_mitogenomes_all.gbff" 2>/dev/null || echo "0")
    echo "  总记录数: ${TOTAL}"
    
else
    echo ""
    echo "⚠️  esearch 未获取到结果（可能是代理问题）"
    echo "   尝试使用 datasets CLI 下载..."
    
    echo ""
    echo "=========================================="
    echo "  备选方案: datasets 下载 Gobiiformes 基因组"
    echo "=========================================="
    datasets download genome taxon "Gobiiformes" \
        --include gbff \
        --filename "${RAW_DIR}/gobiiformes_genomes.zip" \
        2>&1 | tail -5
    
    if [ -f "${RAW_DIR}/gobiiformes_genomes.zip" ]; then
        unzip -o "${RAW_DIR}/gobiiformes_genomes.zip" -d "${RAW_DIR}/ncbi_dataset/" 2>&1 | tail -5
    fi
fi

echo ""
echo "=========================================="
echo "  步骤3: 提取线粒体 CDS 信息概览"
echo "=========================================="
GBFF_FILE="${RAW_DIR}/gobi_mitogenomes_all.gbff"
if [ -f "$GBFF_FILE" ]; then
    echo "  GBFF 文件: ${GBFF_FILE}"
    echo "  总记录数: $(grep -c '^LOCUS' ${GBFF_FILE})"
    echo ""
    echo "  含 CDS 的物种数: $(grep -c '^//' ${GBFF_FILE})"
else
    echo "  ⚠️  无 GBFF 文件，检查下载是否成功"
fi

echo ""
echo "=========================================="
echo "  ✅ 下载阶段完成"
echo "=========================================="
