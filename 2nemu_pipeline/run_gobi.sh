#!/bin/bash
# ============================================================
# run_gobi.sh — 批量运行 NeMu (Gobiiformes)
# 所有基因同时并发，按文件总队列控制并发数
# 并发 10 × 线程 20 = 200 核
# ============================================================

CONDA_BASE="/home/zengjl/miniforge3"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate gobi_mutspec
export JAVA_HOME="${CONDA_BASE}/envs/gobi_mutspec/jdk-17.0.14+7"
export PATH="${JAVA_HOME}/bin:${PATH}"

NEMU=/home/zengjl/gitclone/nemu-pipeline-nf/main_nofilter.nf
INPUT=/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/nemu_input
OUTPUT=/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/nemu_output
TAXDUMP=$HOME/.taxonkit
THREADS=20
MAX_JOBS=20
MAX_TARGETS=50

echo "=========================================="
echo "  批量运行: 所有基因同时并发"
echo "  并发 ${MAX_JOBS}, 线程 ${THREADS}/任务"
echo "  开始: $(date)"
echo "=========================================="

# 收集所有未完成的输入文件（打乱顺序，让各基因同时跑）
ALL_FILES=()
for gene in A6 A8 CO1 CO2 CO3 Cytb ND1 ND2 ND3 ND4 ND4L ND5 ND6; do
    for seq in "${INPUT}/${gene}"/*.fasta; do
        [ -f "$seq" ] || continue
        gene_species=$(basename "$seq" .fasta)
        species="${gene_species#*__}"
        workdir="${OUTPUT}/${gene}/${gene_species}"
    local_work="/mnt/SSD_R0/tmp_zengjl/mtdna-gobiiformes-mutspec/${gene}_${gene_species}"
        [ -f "${workdir}/ms12syn.tsv" ] && continue
        ALL_FILES+=("$seq")
    done
done

# 打乱顺序
RANDOM=$$
for i in $(seq ${#ALL_FILES[@]} -1 1); do
    j=$((RANDOM % i))
    tmp=${ALL_FILES[$i]}
    ALL_FILES[$i]=${ALL_FILES[$j]}
    ALL_FILES[$j]=$tmp
done

TOTAL=${#ALL_FILES[@]}
echo "  待处理: ${TOTAL} 个文件"
echo "=========================================="

COUNT=0
for seq in "${ALL_FILES[@]}"; do
    gene=$(basename "$(dirname "$seq")")
    gene_species=$(basename "$seq" .fasta)
    species="${gene_species#*__}"
    workdir="${OUTPUT}/${gene}/${gene_species}"
    local_work="/mnt/SSD_R0/tmp_zengjl/mtdna-gobiiformes-mutspec/${gene}_${gene_species}"
    
    COUNT=$((COUNT + 1))
    mkdir -p "${workdir}"
    
    echo "  🚀 [${COUNT}/${TOTAL}] ${gene}/${gene_species}"
    
    nextflow -bg -q run "${NEMU}" \
        -w "${local_work}" \
        -output-dir "${workdir}" \
        --input "$(realpath "$seq")" \
        --inputType protein \
        --speciesName "${species}" \
        --gencode 2 \
        --db "/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/ncbi_db/${gene}_BLAST" \
        --taxdump "${TAXDUMP}" \
        --threads ${THREADS} \
        --minSeqs 1 --maxTargetSeqs ${MAX_TARGETS} \
        --model "GTR+FO+G6+I" --modelAsr "GTR+FO+G6+I" \
        --runTreeShrink true --probaArg true \
        --uncertaintyCoef false --consCatCutoff 1 \
        --plot false --spectraType syn \
        > "${workdir}/nemu.log" 2>&1
    
    # 每 MAX_JOBS 个任务暂停等待
    if [ $((COUNT % MAX_JOBS)) -eq 0 ]; then
        echo "  轮次 $((COUNT / MAX_JOBS)), 等待中... ($(date))"
        sleep 120
    fi
    sleep 1
done

echo "等待所有任务完成..."
wait

DONE=$(find "${OUTPUT}" -name "ms12syn.tsv" -not -path "*/work/*" 2>/dev/null | wc -l)
echo "=========================================="
echo "  完成: $(date)"
echo "  成功: ${DONE}/${TOTAL}"
echo "=========================================="
