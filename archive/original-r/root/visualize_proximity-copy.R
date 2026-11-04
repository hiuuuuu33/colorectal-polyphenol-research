# ============================================================
# 网络邻近度分析可视化
# 展示 z_closest 最显著（最负）的化合物
# ============================================================

library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

setwd("/Users/hiugo/中药/数据分析")

# ---- 1. 读取数据 ----
df <- read.csv("merged_zscore_proximity.csv", stringsAsFactors = FALSE)
cat("总记录数:", nrow(df), "\n")
cat("列名:", paste(colnames(df), collapse = ", "), "\n")

# ---- 2. 数据清洗与排序 ----
# 按 z_closest 升序（最负 = 最接近）
df <- df %>%
  mutate(
    chemical = gsub('^"|"$', '', chemical),  # 去引号
    chemical = trimws(chemical)
  ) %>%
  arrange(z_closest)  # 升序：最负的排最前

# Top 10 最接近的化合物
top10 <- head(df, 10)
cat("\n=== Top 10 化合物 (z_closest 最小) ===\n")
top10 %>% select(chemical, z_closest, z_shortest, n_mapped_chemical) %>% print()

# ---- 3. 水平棒棒糖图 (Lollipop) ----
# z_closest 越负 = 越接近 = 越显著
# 颜色按 z_closest 渐变：深色 = 更显著

top10 <- top10 %>%
  mutate(
    chemical = factor(chemical, levels = rev(chemical)),  # 最显著的在上
    neg_z    = -z_closest,                                 # 转为正值便于可视化
    label_text = sprintf("%.2f", z_closest)
  )

p_lollipop <- ggplot(top10, aes(x = neg_z, y = chemical)) +
  # 茎线
  geom_segment(
    aes(xend = 0, yend = chemical),
    linewidth = 1.2, color = "grey60"
  ) +
  # 端点圆
  geom_point(
    aes(size = n_mapped_chemical, fill = neg_z),
    shape = 21, stroke = 1.5, color = "white"
  ) +
  # 直接标注 z 值
  geom_text(
    aes(label = label_text, x = neg_z + max(neg_z) * 0.08),
    size = 3.6, color = "grey30", hjust = 0
  ) +
  # 渐变填充
  scale_fill_gradient(
    low  = "#FEE0D2",    # 浅（不显著）
    high = "#DE2D26",    # 深（显著）
    guide = guide_colorbar(
      title   = expression(-z[closest]),
      reverse = TRUE,
      barwidth  = 0.8,
      barheight = 6
    )
  ) +
  scale_size_continuous(
    range = c(4, 10),
    guide = guide_legend(title = "N targets", override.aes = list(color = "grey50"))
  ) +
  labs(
    title    = "Top 10 Compounds Closest to Carcinogenesis Genes",
    subtitle = "Network proximity z-score (more negative = closer proximity)",
    x        = expression(-z[closest] ~ " (inverted for display)"),
    y        = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    axis.text.y      = element_text(size = 11, face = "italic"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

print(p_lollipop)
ggsave("top10_zclosest_lollipop.png", p_lollipop, width = 12, height = 6, dpi = 300)
cat("\n  top10_zclosest_lollipop.png 已保存\n")

# ---- 4. 水平柱状图 ----
top10 <- top10 %>% mutate(chemical = factor(chemical, levels = rev(levels(chemical))))

p_bar <- ggplot(top10, aes(x = neg_z, y = chemical, fill = neg_z)) +
  geom_col(width = 0.55) +
  geom_text(
    aes(label = label_text, x = neg_z + max(neg_z) * 0.06),
    hjust = 0, size = 3.8, color = "grey20"
  ) +
  scale_fill_gradient(
    low  = "#FEE0D2",
    high = "#DE2D26",
    guide = "none"
  ) +
  labs(
    title    = "Top 10 Compounds Closest to Carcinogenesis Genes",
    subtitle = "Network proximity z-score (more negative = closer proximity)",
    x        = expression(-z[closest] ~ " (inverted for display)"),
    y        = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "grey40", size = 10),
    axis.text.y        = element_text(size = 11, face = "italic"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(p_bar)
ggsave("top10_zclosest_bar.png", p_bar, width = 12, height = 6, dpi = 300)
cat("  top10_zclosest_bar.png 已保存\n")

# ---- 5. 全部化合物的 z_closest 分布 (散点图) ----
p_dist <- ggplot(df, aes(x = reorder(chemical, z_closest), y = z_closest)) +
  geom_point(
    aes(fill = z_closest),
    shape = 21, size = 3, stroke = 0.5, color = "grey70", alpha = 0.85
  ) +
  geom_hline(yintercept = -1.96, linetype = "dashed", color = "#DE2D26", linewidth = 0.5) +
  geom_hline(yintercept =  1.96, linetype = "dashed", color = "#DE2D26", linewidth = 0.5) +
  annotate("text", x = 5, y = -2.2, label = "z = -1.96", size = 3.2, color = "#DE2D26", hjust = 0) +
  annotate("text", x = 5, y =  2.2, label = "z =  1.96", size = 3.2, color = "#DE2D26", hjust = 0) +
  scale_fill_gradient2(
    low  = "#3C5488",
    mid  = "grey90",
    high = "#DC0000",
    midpoint = 0
  ) +
  labs(
    title    = "Network Proximity z-score Distribution",
    subtitle = paste0(nrow(df), " compounds vs carcinogenesis genes\n",
                      "Closest = more negative z-score"),
    x        = "Compound (ordered by z_closest)",
    y        = expression(z[closest])
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position  = "none"
  )

print(p_dist)
ggsave("zscore_distribution.png", p_dist, width = 14, height = 5, dpi = 300)
cat("  zscore_distribution.png 已保存\n")

# ---- 6. Top 20 汇总表格图 (z_closest + targets) ----
top20 <- head(df, 20) %>%
  mutate(rank = row_number()) %>%
  select(rank, chemical, z_closest, z_shortest, n_mapped_chemical)

cat("\n=== Top 20 完整排名表 ===\n")
print(top20, row.names = FALSE)

# 输出 CSV
write.csv(top20, "top20_zclosest.csv", row.names = FALSE)
cat("  top20_zclosest.csv 已保存\n")

cat("\n============================================================\n")
cat("  可视化完成！\n")
cat("============================================================\n")
