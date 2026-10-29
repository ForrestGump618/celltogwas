# celltogwas 1.0.0 升级验收报告

日期：2026-10-29。上游基线：`3f9f36d`。项目：[ForrestGump618/celltogwas](https://github.com/ForrestGump618/celltogwas)。

## 完成内容

原单文件 0.1 工具重构为模块化 Python 包：坐标统一、GTF/GFF3 解析、加权区间索引、表达特异性、调控区域输入、探索性富集、诊断清单和四个 CLI 命令。保留历史 API 与 `gwas_annotaion` 拼写兼容入口。

新增方法包括 Tau/归一化熵/MAD 特异性、距离衰减权重、最大覆盖区间扫描、FFT 环移、染色体流式检验和 BH-FDR。`p_upper` 保留背景，只过滤证据；低绝对表达通过 `min_expression` 过滤。

加入真实 snRNA-seq 衍生 pseudobulk + GWAS 验证、MAGMA/H-MAGMA 实际运算、排名比较和官方批处理/合并脚本。下载清单记录所有公共输入的大小与 SHA-256。

## 本地实测验收

| 检查 | 结果 |
| --- | --- |
| 自动化测试 | 49 项通过；4 个预期的小置换数警告 |
| 核心库单元测试覆盖率 | 83.74% |
| Ruff / 格式 / Git diff 检查 | 通过 |
| wheel 安装 | 全新 Python 3.12 环境，从 site-packages 导入；依赖一致性通过 |
| specificity / validate / annotate / enrich | 均实际执行通过 |
| 安装代码来源 | 安装后核心文件 SHA-256 与真实数据验证源码逐项一致 |
| 真实主分析复现 | 完整 14,161,368 SNP 背景，最终核心重算表逐字节一致 |
| 真实方法比较 | MAGMA v1.10、H-MAGMA 与 celltogwas 均运行；表格、命令和日志随证据包提供 |

核心库覆盖率与验证脚本覆盖率分别报告：脚本的外部运算由真实执行证据验证，不能将真实数据运行伪称为单元测试覆盖率。完整覆盖率表见 `verification.log`。

合成性能基准使用 200,000 SNP、5,000 基因、10 个细胞类与 10,000 次置换。本次机器实测构建约 1.41 秒、注释约 1.00 秒、富集约 0.42 秒；完整参数与运行环境在证据包 benchmark JSON 中。该基准只说明本机工程性能。

## 真实数据结论与边界

完整结果见 `REAL_DATA_VALIDATION.zh-CN.md`。主分析 0/24 类通过 FDR 0.05；全 SNP 敏感性分析 2/24 类通过。未依据结果改变主分析。

只验证了一个 AD accession 与成人 DLPFC。使用公开处理后的真实 pseudobulk，未重跑受控 FASTQ；没有建立跨性状准确率、LD 校准、遗传力分配或因果结论。
外部 MAGMA 实测在 Windows v1.10 完成。仓库包含 Windows/Linux、Python 3.10–3.12 的 CI 定义；远程 CI 状态以 GitHub 当前运行结果为准。

## 发布与复现

发行文件包括 wheel、sdist、源码 ZIP、真实验证证据 ZIP、CLI 演示证据 ZIP、两份中文报告、科学图和 SHA-256 清单。
大型参考数据、MAGMA 可执行文件以及原始受控数据不随发行包分发。MAGMA 按其官方条款从来源下载。

安装 wheel 后运行 `celltogwas --version`。真实验证复现方式见源码 `validation/README.md`；四步 CLI 及安装来源证明保存在演示证据包中。
