# celltogwas 1.0.0 真实数据验证与方法一致性报告

日期：2026-10-29。工程基线：`3f9f36d`。外部验证在 Windows / Python 3.12 上运行。

## 主要结论

完成真实 snRNA-seq 衍生 pseudobulk 与 AD GWAS 的外部示例验证，并实际运行 MAGMA v1.10 和成人脑 H-MAGMA。
预设主分析有 **0/24** 个细胞类通过 BH-FDR 0.05。全 SNP 敏感性分析有 **2/24** 个通过，分别为 Immune (q=0.0023998), Micro.PVM (q=0.031197)。
这些结果支持对免疫相关细胞程序开展后续研究；主分析的统计证据仍有限。结果不足以建立细胞类型的因果作用或跨疾病有效性。

## 数据与质控

| 输入 | 本次实测 |
| --- | --- |
| PsychENCODE / brainSCOPE | 24 个成人 DLPFC 细胞类；公开供者级 log-normalized pseudobulk |
| 供者覆盖 | 并集 352 人；每类 23–346 人；全部 24 类共同供者 1 人 |
| 基因覆盖 | 并集 30,149；全部矩阵均报告 17,158 |
| 标记集合 | 每类 1,000；共 24,000 条，涉及 5,072 个唯一基因 |
| GWAS accession | GCST90027158；欧洲人群；GRCh38 |
| GWAS 原始记录 | 21,101,114 |
| 有效记录中的重复 rsID 行 | 排除 61,103 |
| GRCh37 参考匹配 | 14,161,368 个唯一常染色体 SNP（有效唯一 rsID 的 67.31%） |
| 每 SNP 样本量 | N=n_cases+n_controls；最小 42,140，中位 485,999，最大 487,511 |
| MHC | 主分析包含 78,177 个 chr6:25–35 Mb SNP；另作排除敏感性分析 |

GWAS accession 元数据引用的是 Bellenguez 2020 预印本（10.1101/2020.10.01.20200659），不能将本文件直接标为 2022 最终论文的完整样本。
官方 GWAS MD5 为 `9d23b9ba23532da38ab83fb061bab18f`，本地下载校验一致。1000G EUR 的官方 MD5/SHA1 也一致。全部公共输入的 URL、大小和 SHA-256 记录在 `validation/data_manifest.json`。

原始 FASTQ 涉及受控访问，本次使用开放的真实测序衍生 pseudobulk，未重新处理 FASTQ 或执行细胞级 doublet/QC。没有假定每类都有网站概述的 388 个样本。
仅在所有 24 类均有报告、且可映射到 GENCODE v19 的基因中选择标记，避免将缺行当成零表达。不同细胞类仍存在供者组成及检测率差异。

## 固定分析方案

标记在查看 GWAS 关联结果前固定：Tau 特异性，每类前 1,000，最低最大平均表达 0.1、该类供者检测比例 0.1。输入值保持来源的 log-normalized 尺度，没有反推出原始计数。

GWAS 按唯一 rsID 精确匹配 1000G EUR BIM 的 GRCh37 坐标；不使用静默 liftover。分析使用 P 值，不涉及效应方向或等位基因符号协调。
保留 `N >= 50`，与 MAGMA 默认样本量门槛一致；未匹配和重复记录均有计数。

celltogwas 采用 GENCODE v19 基因体双侧 10 kb、uniform 核、重叠区间最大权重。
主分析对 `P <= 5e-8` 的证据赋 `max(-log10(P)-(-log10(5e-8)), 0)`；阈值只改变证据权重，**完整匹配 SNP 背景保留用于环移**。
执行 10,000 次染色体内独立环移，root seed=13，使用独立子种子，最小偏移比例 0.01；BH 校正在 24 个细胞类内进行。
另运行相同成员的二值集合、`P <= 1e-5`、全 SNP 和排除 MHC 的敏感性分析。

MAGMA 与 H-MAGMA 均由同一个官方 v1.10 可执行文件运行，使用相同 SNP/P/N、EUR 参考和 Ensembl 标记集合。
H-MAGMA 是成人 DLPFC Hi-C SNP→gene 注释，统计引擎仍为 MAGMA；它不是第二个独立统计引擎。
常规 MAGMA 使用 H-MAGMA 作者提供的 `MAGMAdefault.genes.annot`；其基因模型与链方向上下游窗口不同于本次 celltogwas 的 GENCODE ±10 kb，不能当作相同区间的数值等价试验。

## 主分析结果

| 细胞类 | 环移 Z | 经验 P | BH q |
| --- | --- | --- | --- |
| Immune | 2.696 | 0.027597 | 0.54415 |
| Endo__VLMC | 1.755 | 0.062894 | 0.54415 |
| SMC | 1.744 | 0.071993 | 0.54415 |
| PC | 1.485 | 0.090691 | 0.54415 |
| Micro.PVM | 0.990 | 0.12799 | 0.61434 |
| Oligo | 0.010 | 0.32047 | 0.92892 |
| L5.6.NP | -0.168 | 0.39996 | 0.92892 |
| Lamp5 | -0.216 | 0.41586 | 0.92892 |

环移 Z 是相对经验零分布的标准化分数；显著性由经验 P/BH q 判定，不能把 Z 直接代入标准正态分布获得 P。
敏感性结果分别保留，未依据结果更换主分析。各敏感性表的 q 只校正本表 24 类，不代表所有敏感性试验联合控制了错误率。

## 方法一致性

基于共同细胞类，celltogwas 使用环移 Z，MAGMA/H-MAGMA 细胞基因集使用回归 BETA/SE 作排名统计量。

| 比较 | Spearman rho | 描述性 bootstrap 95% 区间 | Top-3 重叠 | Top-5 重叠 |
| --- | --- | --- | --- | --- |
| primary_weighted_vs_standard | 0.5070 | [0.077, 0.791] | 2 | 4 |
| primary_weighted_vs_hmagma | 0.5061 | [0.025, 0.802] | 2 | 4 |
| membership_aligned_vs_standard | 0.5235 | [0.096, 0.806] | 2 | 4 |
| membership_aligned_vs_hmagma | 0.5313 | [0.070, 0.819] | 2 | 4 |
| all_snps_sensitivity_vs_standard | 0.6374 | [0.263, 0.870] | 3 | 4 |
| all_snps_sensitivity_vs_hmagma | 0.5600 | [0.107, 0.848] | 3 | 4 |
| standard_vs_hmagma | 0.8270 | [0.606, 0.933] | 3 | 5 |
| primary_vs_membership_aligned | 0.9670 | [0.867, 0.997] | 3 | 5 |

各方法均在本次 24 个细胞类内单独进行 BH 校正；方法间及敏感性分析间不作联合确认性检验。

| 方法 | BH q<0.05 细胞类 |
| --- | --- |
| celltogwas_primary | 无 |
| celltogwas_membership | 无 |
| celltogwas_all_snps_sensitivity | Immune, Micro.PVM |
| standard_magma | Immune, Micro.PVM |
| hmagma | 无 |

MAGMA/H-MAGMA 共同有 53,021 个含有效 SNP 的基因，基因级 ZSTAT Spearman=0.6009，P 排名 Spearman=0.6009。

Top 基因表按各模型原生基因全集截取；共同基因背景的结果另外保存在 JSON 中。

| Top 基因数 | 重叠数 | Jaccard |
| --- | --- | --- |
| 100 | 40 | 0.2500 |
| 500 | 278 | 0.3850 |
| 1000 | 570 | 0.3986 |

`known_ad_genes_magma.tsv` 列出 APOE、TREM2、BIN1、CR1、PICALM 等预先指定 AD 基因的映射、统计量与共同基因背景中的排名。少数熟知基因的符合情况只作生物学合理性检查。
bootstrap 对基因或细胞类重采样；LD、重叠标记及相关细胞程序不满足独立性，这些区间仅描述当前数据的排名稳定性。不对依赖基因或细胞类的排名相关给出独立样本显著性检验 P 值。

## 可复现性与适用范围

最终核心代码在全部匹配 SNP 上重新计算主分析，结果表逐字节一致：`484580c23e8634be81139a4e054d79780ae60b3c8402f40c3587a31455e8a400`。
发布证据包含完整命令、日志、原始 MAGMA 结果、表格、图及清单。大型公共数据、受控 FASTQ 和 MAGMA 二进制不随包分发。

这里只验证了一个 AD accession 与一种脑区。celltogwas 不直接建模 LD，不估计遗传力，不作 fine-mapping；SNP 索引环移保留排序，但不保证物理碱基距离。
MAGMA/H-MAGMA 也有各自的基因映射与参考人群假设。当前结果适合探索性研究，仍需其他性状、组织、人群及更严格的零假设校准来评估推广性。

复现命令见仓库 `validation/README.md`。主要来源：[brainSCOPE](https://brainscope.gersteinlab.org/)、[pseudobulk 矩阵](https://brainscope.gersteinlab.org/output-pseudobulk-expr.html)、[GWAS accession 原始文件与元数据](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/)、[MAGMA 官方说明](https://cncr.nl/research/magma/)、[H-MAGMA 作者仓库](https://github.com/thewonlab/H-MAGMA)、[GENCODE v19](https://www.gencodegenes.org/human/release_19.html)。
