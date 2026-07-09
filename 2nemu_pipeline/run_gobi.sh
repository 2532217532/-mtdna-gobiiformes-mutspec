#!/bin/bash
# ============================================================
# run_gobi.sh — 批量运行 NeMu (Gobiiformes)
# 并发 10 × 线程 20 = 200 核 (≤218)，降低 NFS IO 争抢
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
MAX_JOBS=10
MAX_TARGETS=50

GENES=("Cytb" "CO1" "CO2" "CO3" "ND1" "ND2" "ND3" "ND4" "ND4L" "ND5" "ND6" "A6" "A8")

echo "=========================================="
echo "  批量运行: 并发 ${MAX_JOBS}, 线程 ${THREADS}/任务"
echo "  开始: $(date)"
echo "=========================================="

PENDING=0

for gene in "${GENES[@]}"; do
    for seq in "${INPUT}/${gene}"/*.fasta; do
        [ -f "$seq" ] || continue
        gene_species=$(basename "$seq" .fasta)
        species="${gene_species#*__}"
        workdir="${OUTPUT}/${gene}/${gene_species}"
        
        [ -f "${workdir}/ms12syn.tsv" ] && continue
        
        PENDING=$((PENDING + 1))
        mkdir -p "${workdir}"
        
        echo "  🚀 #${PENDING} ${gene}/${gene_species}"
        
        nextflow -bg -q run "${NEMU}" \
            -w "${workdir}/work" \
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
        
        # 并发控制: 每启动 MAX_JOBS 个任务，等待一批完成
        if [ $((PENDING % MAX_JOBS)) -eq 0 ]; then
            echo "  等待中... ($(date))"
            sleep 120
        fi
        sleep 1
    done
done

# 等待全部结束
echo "等待所有任务完成..."
wait

echo "=========================================="
echo "  完成: $(date)"
echo "  总任务: ${PENDING}"
echo "  成功: $(find ${OUTPUT} -name ms12syn.tsv 2>/dev/null | wc -l)"
echo "=========================================="
