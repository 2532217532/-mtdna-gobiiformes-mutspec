#!/bin/bash
# 激活 gobi_mutspec 环境
source /home/zengjl/miniforge3/etc/profile.d/conda.sh
conda activate gobi_mutspec

# 设置 JDK17（nextflow 需要）
export JAVA_HOME="/home/zengjl/miniforge3/envs/gobi_mutspec/jdk-17.0.14+7"
export PATH="$JAVA_HOME/bin:$PATH"

echo "✅ gobi_mutspec 环境已激活"
echo "   Python:  $(python --version 2>&1)"
echo "   Nextflow: $(nextflow -version 2>&1 | head -2 | tail -1)"
echo "   datasets: $(datasets version 2>/dev/null | head -1)"
echo "   BLAST:    $(blastn -version 2>&1 | head -1)"
