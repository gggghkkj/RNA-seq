# 复现图4：形态多能干细胞转录组分析

## 项目概述

本项目复现论文*"Formative pluripotent stem cells show features of epiblast cells poised for gastrulation"*（Cell Research, 2021）中的核心转录组分析（图4a-c）。使用来自小鼠胚胎干细胞（mESCs）、上胚层干细胞（EpiSCs）和形态多能干细胞（fPSCs）的RNA-seq数据，进行差异表达分析、PCA和可视化，以表征这三种多能性状态的转录组特征。

## 分析流程

```
原始数据 (FASTQ) 
    ↓
质量控制 (FastQC + MultiQC)
    ↓
接头修剪 (fastp)
    ↓
参考基因组索引构建 (HISAT2)
    ↓
测序reads比对 (HISAT2 → samtools)
    ↓
基因定量 (featureCounts)
    ↓
差异表达分析 (DESeq2)
    ↓
可视化 (PCA, 热图, 点图)
```

分步说明：

1. 质量控制 — FastQC评估原始测序reads质量；MultiQC生成汇总报告
2. 接头修剪 — fastp去除接头并过滤低质量reads（Q20，长度≥50bp）
3. 基因组索引构建 — HISAT2为GRCm39（mm10）构建基因组索引
4. 测序reads比对 — HISAT2将双端reads比对到小鼠基因组（比对率≥95%）
5. 基因定量 — featureCounts使用Ensembl GTF注释计算每个基因的reads数
6. 差异表达分析 — DESeq2进行标准化（VST）、成对比较并鉴定差异表达基因（DEGs）
7. 可视化 — PCA图（图4a）、前50个高变异基因热图（图4b）、标记基因表达点图（图4c）

## 环境与依赖

| 软件 | 版本 | 用途 |
|------|------|------|
| HISAT2 | 2.2.3 | 测序reads比对 |
| samtools | 1.24 | BAM文件处理 |
| featureCounts | 2.1.1 | 基因定量 |
| fastp | 1.3.6 | 接头修剪 |
| FastQC | 0.12.1 | 质量控制 |
| MultiQC | 1.35 | QC报告汇总 |
| R | 4.5.3 | 统计分析 |
| DESeq2 | 1.50.2 | 差异表达分析 |
| ggplot2 | 4.0.3 | 可视化 |
| pheatmap | 1.0.13 | 热图生成 |

Conda环境设置：

```bash
# 创建conda环境
mamba create -n rnaseq -c bioconda -c conda-forge \
    hisat2 samtools subread fastp fastqc multiqc r-base

# 安装R包
mamba install -n rnaseq -c bioconda -c conda-forge \
    bioconductor-deseq2 bioconductor-annotationdbi \
    r-ggplot2 r-pheatmap r-rcolorbrewer r-circlize

# 激活环境
conda activate rnaseq
```

## 运行方法

### 前提条件

1. 下载原始数据（fastq.gz）到`data/raw/`
2. 下载参考基因组和注释文件到`data/raw/`：
   - `Mus_musculus.GRCm39.dna.primary_assembly.fa.gz` 来自 [Ensembl](https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/)
   - `Mus_musculus.GRCm39.113.gtf.gz` 来自 [Ensembl](https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/)
3. 解压参考文件：
   ```bash
   gunzip data/raw/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
   gunzip data/raw/Mus_musculus.GRCm39.113.gtf.gz
   ```

### 运行完整流程

```bash
# 执行完整分析流程
bash run_pipeline.sh
```

### 运行单个步骤

```bash
# 步骤1-5：Linux命令（质量控制、修剪、比对、定量）
# 这些命令嵌入在run_pipeline.sh中

# 步骤6：DESeq2分析（R）
Rscript scripts/01_deseq2_analysis.R

# 步骤7：生成图表（R）
Rscript scripts/02_figures.R
```

### 输出文件

| 文件 | 描述 |
|------|------|
| `results/figures/Figure4a_PCA.pdf` | 显示样本聚类的PCA图 |
| `results/figures/Figure4b_Heatmap.pdf` | 前50个高变异基因热图 |
| `results/figures/Figure4c_Markers_dotplot.pdf` | 标记基因表达点图 |
| `results/tables/DEGs_*.csv` | 差异表达基因列表 |
| `data/processed/FPKM_matrix.csv` | 标准化FPKM表达矩阵 |

## 数据来源

RNA-seq数据下载自GEO Series GSE154290：

[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154290)

本分析中使用的样本：

| 样本 | GEO登录号 | SRA登录号 | 细胞类型 |
|------|-----------|-----------|----------|
| mESCs rep1 | GSM4668511 | SRR12207279 | 原始态mESCs |
| mESCs rep2 | GSM4668512 | SRR12207280 | 原始态mESCs |
| EpiSCs rep1 | GSM4668515 | SRR12207283 | 始发态EpiSCs |
| EpiSCs rep2 | GSM4668516 | SRR12207284 | 始发态EpiSCs |
| fPSCs P1 | GSM4668520 | SRR12207288 | 形态多能干细胞（第1代） |
| fPSCs P10 | GSM4668521 | SRR12207289 | 形态多能干细胞（第10代） |
| fPSCs P20 | GSM4668522 | SRR12207290 | 形态多能干细胞（第20代） |
| fPSCs P30 | GSM4668523 | SRR12207291 | 形态多能干细胞（第30代） |

下载示例：

```bash
# 使用SRA Toolkit
prefetch SRR12207279
fastq-dump --split-files SRR12207279

# 或使用wget（双端测序）
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12207279/SRR12207279_1.fastq.gz
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12207279/SRR12207279_2.fastq.gz
```

## 参考基因组

- 组装版本：GRCm39（mm10）— Mus musculus
- Ensembl版本：113
- 下载：[Ensembl FTP](https://ftp.ensembl.org/pub/release-113/)

## 引用

如果您使用此代码或分析流程，请引用：

> Wang X, Xiang Y, Yu Y, et al. Formative pluripotent stem cells show features of epiblast cells poised for gastrulation. *Cell Research*. 2021;31:526-541. doi:[10.1038/s41422-021-00477-x](https://doi.org/10.1038/s41422-021-00477-x)

## 许可证

本项目用于教育和研究目的。
