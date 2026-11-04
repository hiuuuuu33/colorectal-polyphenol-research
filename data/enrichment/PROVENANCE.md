# Frozen enrichment results

The 15 CSV files under `frozen/` are original, complete GO/KEGG result-table exports for the three contrasts. They retain the term descriptions, gene ratios, background ratios, gene lists, raw P values, and adjusted P values used for the submitted analysis. No statistical values are replaced during figure generation.

Original locations within the supplied project folder:

- Adenoma versus normal: `数据分析/normal_vs_adenoma/`.
- CRC versus adenoma: `数据分析/adenoma_vs_cancer/`.
- CRC versus normal: `数据分析/` (the five contrast-specific CSV files).

`scripts/R/02_enrichment_figures.R` selects rows with `p.adjust < 0.05`, exports the significant tables, verifies all reported term counts, and redraws the top 15 terms per panel. Figure 4 KEGG panels therefore include all 14 upregulated pathways and the single downregulated pathway. The full tables remain available for every contrast.

This is a frozen-result reproduction. It regenerates figures and summaries from the saved original enrichment output without contacting a changing online annotation service. The original analysis used clusterProfiler defaults for the enrichment background and separate up/down KEGG analyses. A fresh online enrichment calculation may return different results as annotations change; the original CSVs are retained to keep the submitted figures traceable.
