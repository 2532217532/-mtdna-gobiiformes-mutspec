# Gobiiformes Mitochondrial Mutational Spectra

使用 NeMu 管道从 NCBI 下载的 Gobiiformes（虾虎鱼目）mtDNA 数据推导 12/192 组分突变谱，
并进行 COSMIC 签名分配分析。

## 项目结构

```
mtdna-gobiiformes-mutspec/
├── 1data_preparation/       ← 数据下载、CDS 提取、NeMu 输入准备
│   ├── 0download_ncbi.sh    # NCBI datasets 下载脚本（待编写）
│   ├── 1extract_cds.ipynb   # 从 GenBank 提取 13 个 mt 基因 CDS（待编写）
│   ├── 2prepare_input.ipynb # 生成 NeMu 输入格式 + BLAST DB（待编写）
│   └── markdb.sh            # BLAST 数据库构建脚本（已复制）
│
├── 2nemu_pipeline/          ← NeMu 突变谱计算管道
│   ├── nemu.nf              # NeMu 核心 Nextflow 管道
│   ├── nemu_goby.config     # Gobiiformes 管道参数配置
│   ├── batch_nemu.nf        # 批量运行 Nextflow 脚本
│   └── run_gobi.sh          # bash 批量运行脚本
│
├── 3collect_spectra/        ← 聚合 NeMu 输出为最终数据集
│   └── collect_spectra.ipynb
│
├── 4signatures/             ← COSMIC 签名分配
│   ├── 1signatures_analysis_sigpro.ipynb  # SigProfilerAssignment
│   ├── 3mSigAct_analysis.R                # mSigAct 签名分解
│   └── ...
│
├── data/
│   ├── raw/                 # NCBI 下载原始数据（gitignored）
│   ├── nemu_input/          # NeMu 输入 FASTA（gitignored）
│   ├── nemu_output/         # NeMu 输出结果（gitignored）
│   └── dataset/             # 最终聚合数据集
│
├── requirements.txt         # Python 依赖
├── .gitignore
└── README.md
```

## 流程概览

1. **数据准备**: NCBI datasets 下载 Gobiiformes 基因组 → 提取 13 mt 基因 CDS → 准备蛋白序列 + BLAST DB
2. **NeMu 管道**: 每个 gene-species 独立运行，得到 12/192 组分突变谱
3. **结果聚合**: 汇总所有 gene-species 的突变谱为最终数据集
4. **签名分配**: SigProfilerAssignment + mSigAct 分解 COSMIC 签名

## 依赖

- Python 3.9+
- Nextflow
- IQ-TREE2
- MAFFT / MACSE
- NCBI BLAST+
- NeMu-pipeline（Singularity 容器）
- R 4.4.1（mSigAct）
- NCBI datasets CLI

详见 `requirements.txt`。

## 数据来源

NCBI GenBank，通过 `datasets` CLI 下载 Gobiiformes（taxid: 1489878）的参考基因组。
