#!/bin/bash
# 批量构建 BLAST V5 数据库（适用于 Gobiiformes NCBI 数据）
# 用法: 先确保 data/ncbi_db/ 目录下每个基因有对应的 nt FASTA

mkdir -p /home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/ncbi_db

# 线粒体基因列表
GENES=("Cytb" "CO1" "CO2" "CO3" "ND1" "ND2" "ND3" "ND4" "ND4L" "ND5" "ND6" "A6" "A8")

for gene in "${GENES[@]}"; do
    echo "====================================="
    echo "构建 BLAST V5 库：${gene}"
    
    INPUT_FA="/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/ncbi_db/${gene}.fasta"
    OUT_DB="/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/data/ncbi_db/${gene}_BLAST"

    if [ ! -f "$INPUT_FA" ]; then
        echo "⚠️  文件不存在：${INPUT_FA}，跳过"
        continue
    fi

    makeblastdb -in "${INPUT_FA}" \
                -dbtype nucl \
                -out "${OUT_DB}" \
                -parse_seqids

    echo "${gene} V5 数据库构建 ✅"
done

echo "====================================="
echo "🎉 所有基因库构建完成！"
