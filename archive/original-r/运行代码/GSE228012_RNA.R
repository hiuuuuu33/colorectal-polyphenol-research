# 加载必要的包
# 注意：FPKM数据使用limma分析，不需要DESeq2/edgeR（它们需要原始计数数据）
# library(DESeq2)
# library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)


# 下载补充文件中的FPKM数据，然后转换为计数数据
setwd("d:/Bolejiaoyu/maixiao")

# 设置下载方式：Windows 上 SSL 吊销检查可能失败，使用 curl R 包直连
options(timeout = max(600, getOption("timeout")))

url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228012/suppl/GSE228012_all.gene.FPKM.mRNA_processed.txt.gz"

# 方法1：使用 method="curl" + --ssl-no-revoke 绕过吊销检查
cat("正在下载数据...\n")
download_success <- FALSE

# 先尝试 curl 方法 + --ssl-no-revoke（绕过 SSL 吊销检查）
tryCatch({
  download.file(url, destfile = "GSE228012_FPKM.txt.gz",
                method = "curl", mode = "wb",
                extra = c("--ssl-no-revoke", "-L", "--connect-timeout", "30"))
  download_success <- TRUE
  cat("下载成功 (curl --ssl-no-revoke)\n")
}, error = function(e) {
  cat("curl --ssl-no-revoke 失败:", conditionMessage(e), "\n")
})

# 回退：使用 R 内置的 curl 包
if (!download_success) {
  tryCatch({
    if (!requireNamespace("curl", quietly = TRUE)) install.packages("curl")
    h <- curl::new_handle()
    curl::handle_setopt(h, ssl_verifypeer = FALSE, ssl_verifyhost = FALSE,
                        followlocation = TRUE, connecttimeout = 30)
    curl::curl_download(url, "GSE228012_FPKM.txt.gz", handle = h)
    download_success <- TRUE
    cat("下载成功 (curl R package, SSL verify off)\n")
  }, error = function(e) {
    cat("curl R package 失败:", conditionMessage(e), "\n")
  })
}

# 最后尝试：httr 包
if (!download_success) {
  tryCatch({
    if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
    resp <- httr::GET(url, httr::config(ssl_verifypeer = FALSE),
                       httr::write_disk("GSE228012_FPKM.txt.gz", overwrite = TRUE),
                       httr::timeout(60))
    if (httr::status_code(resp) == 200) {
      download_success <- TRUE
      cat("下载成功 (httr, SSL verify off)\n")
    } else {
      cat("httr HTTP状态:", httr::status_code(resp), "\n")
    }
  }, error = function(e) {
    cat("httr 失败:", conditionMessage(e), "\n")
  })
}

if (!download_success) {
  stop("所有下载方法均失败，请检查网络连接或手动下载文件放置于工作目录")
}
# gunzip("GSE228012_FPKM.txt.gz")
# 直接读取FPKM数据
fpkm_data <- read.table(gzfile("GSE228012_FPKM.txt.gz"),
                        header = TRUE, sep = "\t")



# 查看数据维度
print(dim(fpkm_data))
head(fpkm_data[, 1:5])

# 实际数据格式：第1列 GeneID（数字），第2列 X（基因名），第3列起为FPKM表达值
# 删除 GeneID 列，使数据结构与后续代码一致
fpkm_data$GeneID <- NULL

# 方法1：对于重复基因，保留表达量最高的行

# 查看重复的基因名
dup_genes <- fpkm_data$X[duplicated(fpkm_data$X)]
cat("重复基因数量:", length(dup_genes), "\n")
print(head(dup_genes))

# 按行平均值排序，保留每个基因的最大表达值
fpkm_data$mean_expression <- rowMeans(fpkm_data[, -1], na.rm = TRUE)
fpkm_data <- fpkm_data[order(fpkm_data$mean_expression, decreasing = TRUE), ]

# 去除重复，保留第一个（即表达量最高的）
fpkm_data <- fpkm_data[!duplicated(fpkm_data$X), ]

# 删除辅助列
fpkm_data$mean_expression <- NULL

# 设置行名
rownames(fpkm_data) <- fpkm_data$X
fpkm_data <- fpkm_data[, -1]

cat("去重后基因数量:", nrow(fpkm_data), "\n")



#对于FPKM数据，应该使用limma而不是DESeq2/edgeR，因为limma可以处理连续的表达值，而DESeq2/edgeR需要原始计数数据。
# FPKM转TPM
# 创建样本分组
# 方法2：直接清理列名，然后分组
clean_sample_names <- gsub("\\.FPKM$", "", colnames(fpkm_data))
clean_sample_names <- trimws(clean_sample_names)

# 更新数据框列名
colnames(fpkm_data) <- clean_sample_names

# 现在使用清理后的列名进行分组
sample_names <- colnames(fpkm_data)
groups <- ifelse(grepl("T$", sample_names), "CRC",
                 ifelse(grepl("N$", sample_names), "Normal",
                        ifelse(grepl("pre$", sample_names), "Adenoma", "Unknown")))

sample_info <- data.frame(
  sample = sample_names,
  group = factor(groups, levels = c("Normal", "Adenoma", "CRC")),
  row.names = sample_names
)

print("清理后的样本分组情况:")
print(table(sample_info$group))

# 过滤低表达基因（在至少20%样本中FPKM > 1）
keep_genes <- rowSums(fpkm_data > 1) >= round(0.2 * ncol(fpkm_data))
filtered_fpkm <- fpkm_data[keep_genes, ]
cat("过滤后基因数:", nrow(filtered_fpkm), "\n")

# log2转换
log2_fpkm <- log2(filtered_fpkm + 1)

# limma差异分析
design <- model.matrix(~0 + sample_info$group)
colnames(design) <- c("Normal", "Adenoma", "CRC")

fit <- lmFit(log2_fpkm, design)

cont.matrix <- makeContrasts(
  CRC_vs_Normal = CRC - Normal,
  Adenoma_vs_Normal = Adenoma - Normal,
  CRC_vs_Adenoma = CRC - Adenoma,
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 提取结果
get_de_results <- function(fit, contrast, fdr_threshold = 0.05, logfc_threshold = 1) {
  results <- topTable(fit, coef = contrast, number = Inf, adjust.method = "fdr")
  results$gene <- rownames(results)
  
  sig_genes <- results[results$adj.P.Val < fdr_threshold & 
                         abs(results$logFC) > logfc_threshold, ]
  
  cat(contrast, ":\n")
  cat("  总显著基因:", nrow(sig_genes), "\n")
  cat("  上调基因:", sum(sig_genes$logFC > 0), "\n")
  cat("  下调基因:", sum(sig_genes$logFC < 0), "\n\n")
  
  return(list(full = results, significant = sig_genes))
}

# 获取各比较组结果
crc_normal <- get_de_results(fit2, "CRC_vs_Normal")
adenoma_normal <- get_de_results(fit2, "Adenoma_vs_Normal")
crc_adenoma <- get_de_results(fit2, "CRC_vs_Adenoma")

# 保存中间结果，供下游分析脚本使用
save(log2_fpkm, sample_info, fit2,
     crc_normal, adenoma_normal, crc_adenoma,
     file = "GSE228012_limma_results.Rdata")
cat("中间结果已保存至 GSE228012_limma_results.Rdata\n")
