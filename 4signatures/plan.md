# Plan: 阶段4 — Gobiiformes mtDNA COSMIC 签名分解

## 背景

阶段 1–3 已完成：
- 295 个物种，13 个线粒体基因，3,543 个 NeMu 谱成功
- 输出：`data/processed/species_spectra_192syn.csv`（295 行 × 193 列）及基因级别文件

`4signatures/` 目录代码从 chordata 模板项目 `/home/zengjl/WorkSpace/mtdna-192component-mutspec-chordata/` 复制而来，引用的是旧的脊椎动物数据集，且硬编码了脊椎动物纲名（Actinopteri, Amphibia 等）。需要适配到新的 Gobiiformes NeMu 数据，同时复用已有工具函数。

**关键缺口：** 处理后的谱数据没有 family（科）分类列，需从 `data/raw/gbff/` 的 GenBank GBFF 文件中提取。

**已有资源：** chordata 模板项目已生成 `triplet_counts_GRCh37.json`（7.7 KB），可直接复制，无需重新下载 3 GB 人类基因组并重新计数。路径：
`/home/zengjl/WorkSpace/mtdna-192component-mutspec-chordata/4signatures/data/triplet_counts_GRCh37.json`

**重要差异——链标记方向：**
- 旧 chordata notebook 输入数据使用 L-strand 标记法，notebook 中调用 `rev_comp()` 转为 H-strand
- 我们新的 NeMu 数据**已经是 H-strand**（`A[C>A]T` 格式），加载时**必须跳过 rev_comp 步骤**

**`4signatures/data/` 初始状态：** 完全为空（仅 `.gitignore`）。

---

## 设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 分析粒度 | **物种级别**（13 基因求平均） | 每个物种数据最大化；295 个物种提供稳健的科级均值 |
| 分组方式 | **科（Family）**（从 GBFF 提取） | 与 chordata 模板的纲级别分组逻辑一致；科是虾虎鱼目内自然的分类单元 |
| 数据来源 | **新 NeMu 输出**（`species_spectra_192syn.csv`） | 比旧的 4 基因数据集更丰富；使用全部 13 个 mtDNA 基因 |
| 链标记 | **跳过 `rev_comp()`** | NeMu 输出已是 H-strand；原 notebook 的 `rev_comp()` 是为 L-strand 输入服务的 |
| 签名工具 | **SigProfilerAssignment + mSigAct 双轨** | 交叉验证；与 chordata 论文方法一致 |
| 三联体计数 | **从 chordata 项目复制** | 避免重新下载 3 GB 人类基因组；使用相同的 GRCh37 参考 |
| 科过滤阈值 | **每科 ≥3 个物种** | 确保均值谱稳健；仅 1-2 个物种的小科噪声过大 |

---

## 步骤 1：创建目录结构并复制三联体计数

**目的：** SigProfilerAssignment 需要将非人物种的谱重归一化到人类基因组三核苷酸上下文。计数文件在 chordata 项目中已存在，直接复用。

**操作：**
1. 在 `4signatures/data/` 下创建目录结构：
   ```
   SigProfilerAssignment/input/
   SigProfilerAssignment/output/
   mSigAct/input/
   mSigAct/output/raw_output/
   mSigAct/output/figures/
   ```
2. 从 chordata 项目复制三联体计数：
   ```bash
   cp /home/zengjl/WorkSpace/mtdna-192component-mutspec-chordata/4signatures/data/triplet_counts_GRCh37.json \
      /home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec/4signatures/data/
   ```

**涉及文件：** 目录创建 + 复制已有 JSON（不需要运行 `0count_human_triplets_freqs.py`）。

---

## 步骤 2：从 GBFF 文件提取科级分类信息

**目的：** 物种谱缺少 `family` 列。需要科级分类将物种分组，类似于脊椎动物分析中按纲分组。

**操作：**
创建新脚本 `4signatures/extract_taxonomy.py`：
1. 解析 `data/raw/gbff/` 中全部 7 个批次 GBFF 文件
2. 对每条记录提取：
   - `accession`（VERSION 行）
   - `species`（SOURCE 行中的 organism 名称）
   - 完整谱系（ORGANISM 块中分号分隔的行）
3. 从谱系中提取 **科（family）**（以 `-idae` 结尾的分类阶元）。虾虎鱼目中预计出现的科包括：Gobiidae, Eleotridae, Oxudercidae, Odontobutidae, Butidae, Milyeringidae, Rhyacichthyidae, Thalasseleotrididae 等
4. 输出 `data/processed/taxonomy.csv`，列名为：`species, accession, family, lineage`

**验证：** 打印各科物种计数——确认大部分物种归属 Gobiidae（最大的科）。

**新文件：** `4signatures/extract_taxonomy.py`

---

## 步骤 3：适配并运行 SigProfilerAssignment notebook

**目的：** 核心分析——将 Gobiiformes 突变谱分解为 COSMIC SBS 签名。

**操作：** 创建 `4signatures/1signatures_analysis_gobi.ipynb`（基于 `1signatures_analysis_sigpro.ipynb` 改写）：

### 数据加载与准备
1. 加载 `data/processed/species_spectra_192syn.csv`（宽表格式，295 物种 × 192 SBS 频率）
2. 加载 `data/processed/taxonomy.csv`，按 `species` 列合并，添加 `family` 列
3. **跳过 `rev_comp()`**——NeMu 输出已是 H-strand 标记法（`A[C>A]T` 格式，匹配 `possible_sbs192`）。原 notebook 因输入是 L-strand 才需要 `rev_comp()`
4. 使用 `utils.py` 的 `complete_sbs_columns()` 补全缺失列并重排为标准顺序

### 科级聚合
5. 按 `family` 分组，计算每科**均值**谱（与 `calc_mutspec_class` 逻辑一致）
6. 可选：额外创建 "Gobiiformes_all" 行作为全目级别的谱
7. 过滤掉物种数 <3 的科，确保均值谱稳健

### Low/High/Diff 拆分与重归一化
8. **复用原 notebook cell 13–20 的逻辑：**
   - 将各科谱拆分为 Low-Ts（G>A, T>C）、High-Ts（C>T, A>G）和 Tv（所有颠换，各取半）
   - 创建 `__Ts only` 和 `__Ts & Tv` 两种变体
   - 使用 `human_counts` 和 `save_wide_cls_spectra()` 重归一化到人类基因组
9. 写入输入文件：
   - `4signatures/data/SigProfilerAssignment/input/low_Ts_samples.txt`
   - `4signatures/data/SigProfilerAssignment/input/high_Ts_samples.txt`
   - `4signatures/data/SigProfilerAssignment/input/high_minus_low_Ts_samples.txt`

### SigProfilerAssignment cosmic_fit
10. **复用原参数**（cosmic_version=3.3, 相同的排除签名子组, 相同的 nnls 惩罚参数）。对 Low、High、Diff 各运行 1 次，共 3 次

### 绘图
11. 使用 `plotActivity.py`（原样复用）生成活动堆积柱状图
12. 将绘图 cell 中的科名引用从脊椎动物纲名更新为虾虎鱼目科名

**新文件：** `4signatures/1signatures_analysis_gobi.ipynb`  
**复用：** `utils.py`, `plotActivity.py`（无需修改）

---

## 步骤 4：准备 mSigAct 先验

**目的：** 管道 B 轨使用 mSigAct 进行交叉验证。先验从 SigProfilerAssignment 输出计算。

**操作：** 运行 `4signatures/2prepare_priors_for_mSigAct.py`，仅需更新输入路径指向新的 SigProfiler 输出。该脚本读取 `Assignment_Solution_Activities.txt`，对 SBS 计数求和，计算占比 ≥1% 的签名。输出到 stdout——我们捕获 6–8 个 top SBS 签名及其比例。

**修改文件：** `4signatures/2prepare_priors_for_mSigAct.py`（仅改路径）

---

## 步骤 5：适配并运行 mSigAct 分析（R）

**目的：** mSigAct 提供独立的分解方法进行交叉验证。

**操作：** 创建 `4signatures/3mSigAct_analysis_gobi.R`（基于 `3mSigAct_analysis.R` 改写）：
1. 更新输入路径指向新的 Gobiiformes SigProfilerAssignment 输入文件
2. 用步骤 4 得出的 SBS 签名和比例更新 `sig_use` 和 `sig_prop`
3. 运行 3 种模式：`custom_prop`（SigProfiler 先验）、`prop1`（所有相关 SBS 均匀先验）、`custom_from1`（从 prop1 结果推算先验）
4. 更新输出目录路径

**新文件：** `4signatures/3mSigAct_analysis_gobi.R`

---

## 步骤 6：汇总 mSigAct 输出并绘图

**目的：** 将逐样本的 mSigAct 结果合并为表格，生成论文级别的图。

**操作：**
1. 创建 `4signatures/4aggregate_mSigAct_outputs_gobi.ipynb`（基于原版改写）：
   - 从 `prop1` 运行结果计算 top SBS 比例 → 指导 `custom_from1` 先验
   - 将所有 exposure CSV 透视到宽表活动矩阵
   - 收集所有距离指标
2. 创建 `4signatures/5plot_output_gobi.py`（基于 `5plot_output_for_mSigAct.py` 改写）：
   - 将硬编码的 `set_order`（脊椎动物纲名）替换为虾虎鱼目科名
   - 生成 4 张 PDF 图

**新文件：** `4signatures/4aggregate_mSigAct_outputs_gobi.ipynb`, `4signatures/5plot_output_gobi.py`  
**复用：** `plotActivity.py`（无需修改）

---

## 文件操作汇总

| 操作 | 文件 | 说明 |
|------|------|------|
| 复制 | `4signatures/data/triplet_counts_GRCh37.json` | 从 chordata 项目复制（免去重新生成） |
| 新建 | `4signatures/extract_taxonomy.py` | 从 GBFF 提取科级分类 |
| 新建 | `4signatures/1signatures_analysis_gobi.ipynb` | 主 SigProfilerAssignment 分析 |
| 修改 | `4signatures/2prepare_priors_for_mSigAct.py` | 更新路径指向新 SigProfiler 输出 |
| 新建 | `4signatures/3mSigAct_analysis_gobi.R` | mSigAct 使用虾虎鱼科级先验 |
| 新建 | `4signatures/4aggregate_mSigAct_outputs_gobi.ipynb` | 汇总 mSigAct 结果 |
| 新建 | `4signatures/5plot_output_gobi.py` | 用科名绘制 mSigAct 结果 |
| 复用 | `4signatures/utils.py` | 无需修改 |
| 复用 | `4signatures/plotActivity.py` | 无需修改 |

---

## 验证

1. **分类提取**：打印各科物种计数，手动抽查 5–10 个物种与 NCBI Taxonomy 核对
2. **SigProfilerAssignment**：检查 `Solution_Stats.txt`——余弦相似度应 >0.8（多数科）
3. **mSigAct**：检查距离 CSV——余弦距离应 <0.1（良好拟合）
4. **交叉验证**：对比 SigProfilerAssignment 与 mSigAct 的 top SBS 签名——主导签名应有 ≥80% 重叠
5. **图表**：目视检查堆积柱状图——科间差异应可见

---

## 执行顺序（依赖关系）

```
步骤1（复制三联体计数+创建目录）──┐
                                  ├──> 步骤3（SigProfiler notebook）──> 步骤4（先验）──┐
步骤2（从GBFF提取科级分类）──────┘                                                  │
                                                                                       ├──> 步骤5（mSigAct R）
                                                                                       │
                                                                                       └──> 步骤6（汇总+绘图）
```

步骤 1 和 2 相互独立，可并行进行。  
步骤 3 依赖步骤 1+2。  
步骤 4–6 构成 mSigAct 轨道，依赖步骤 3 的 SigProfilerAssignment 输出。

注意：原有的 `0count_human_triplets_freqs.py` 和 `1signatures_analysis_sigpro.ipynb` 保留不动作为参考。所有新工作写入带 `_gobi` 后缀的新文件。
