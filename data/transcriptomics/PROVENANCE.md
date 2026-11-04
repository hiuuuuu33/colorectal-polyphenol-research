# Transcriptomic inputs

`GSE228012_all.gene.FPKM.mRNA_processed.txt` is a byte-for-byte copy of the processed expression matrix supplied in the project folder. The public accession is [GSE228012](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228012). It contains 18,116 rows and 33 tissue samples. The source file uses Entrez gene identifiers and gene symbols; numeric sample identifiers are preserved as supplied.

The four files under `reference/` are the original exported significant-gene tables and progression-associated gene table. They are read only for verification, not as inputs to the new differential-expression calculation.

The reproduction follows the submitted analysis: retain the highest mean-expression row for duplicate symbols; keep genes with FPKM > 1 in at least seven samples; apply log2(FPKM + 1); fit an unblocked limma group model; apply BH-adjusted P < 0.05 and |log2 fold change| > 1. The two adjacent comparisons define the same 29 progression-associated genes (8 up, 21 down). This is the analysis represented in the submitted paper. PCA uses unit-variance gene scaling to match the original FactoMineR default; eigenvector signs and plot layout can vary across implementations.

Source locations within the original project folder:

| Repository file | Original relative path |
|---|---|
| expression matrix | `GSE228012_all.gene.FPKM.mRNA_processed.txt` |
| `reference/adenoma_normal_diff.csv` | `数据分析/差异分析对比/adenoma_normal_diff.csv` |
| `reference/crc_adenoma_diff.csv` | `数据分析/差异分析对比/crc_adenoma_diff.csv` |
| `reference/crc_normal_diff.csv` | `数据分析/差异分析对比/crc_normal_diff.csv` |
| `reference/progressive_genes.csv` | `数据分析/progressive_analysis/progressive_genes.csv` |

The scripts generate new figures from these inputs. Original illustrations are archived separately; the regenerated figures preserve the data and panel content while using a consistent, portable rendering workflow.
