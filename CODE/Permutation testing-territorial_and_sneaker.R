library(dplyr) 
library(ggplot2) # 用于绘制直方图

## -----------------------------
## Den-based permutation test
## Den-level delta (average of within-Den differences)
## Territorial vs Sneaker
## -----------------------------

## ========== 1) 读入数据 ========== 
## 假设你的数据已经在数据框 df 中，列名包括：
## MaleID, Stream, Season, Den, Tactic, Total length (m),
## Body weight (kg), Body mass condition

df <- read.csv("F:/for_R/paper_mating_tacitcs_for_genetic_diversity/re_ana/tactics_reana/8.25/body_size_of_males.csv",
               header = TRUE,
               na.strings = c("U", "", "NA", "NaN", "."), # 保持与 GLMM 脚本一致的 NA 处理
               check.names = FALSE, # 保持原始列名，以便后续 rename
               stringsAsFactors = FALSE,
               fileEncoding = "UTF-8-BOM")


df <- df %>%
  rename(
    TotalLength = `Total length (m)`,
    BodyWeight = `Body weight (kg)`,
    BodyCondition = `Body mass condition`
  ) %>%
  mutate(
    MaleID = factor(trimws(MaleID)),
    Stream = factor(trimws(Stream)),
    Season = factor(Season),
    Den = factor(trimws(Den)),
    Tactic = factor(
      trimws(Tactic),
      levels = c("Territorial", "Sneaker")
    ),
    Tactic_bin = ifelse(Tactic == "Sneaker", 1L, 0L),
    DenID = factor(
      interaction(Stream, Season, Den, sep = "_", drop = TRUE)
    )
  )


n_perm <- 9999
seed <- 123
set.seed(seed)

resp_list <- c("TotalLength",
               "BodyWeight",
               "BodyCondition")

group_levels <- c("Territorial", "Sneaker")


df2 <- df


df2 <- df2[df2$Tactic %in% group_levels, , drop = FALSE]

df2$Tactic <- as.character(df2$Tactic)
df2$Den    <- as.character(df2$Den)

if (length(unique(df2$Tactic)) < 2) {
  stop("在删除 U 之后，数据中至少需要同时存在 Territorial 和 Sneaker。")
}


den_delta_mean_stat <- function(data, response, den_var = "Den",
                                group_var = "Tactic",
                                g1 = "Territorial", g2 = "Sneaker") {
  
  dens <- unique(data[[den_var]])
  deltas <- numeric(0)
  
  for (d in dens) {
    idx <- data[[den_var]] == d
    y_d <- data[[response]][idx]
    t_d <- data[[group_var]][idx]
    
    if (any(t_d == g1) && any(t_d == g2)) {
      m1 <- mean(y_d[t_d == g1], na.rm = TRUE)
      m2 <- mean(y_d[t_d == g2], na.rm = TRUE)
      deltas <- c(deltas, m1 - m2)
    }
  }
  
  
  mean(deltas, na.rm = TRUE) 
}

## ==========置换检验主函数：Den 内置换 tactic ========== 
perm_test_by_den_delta_mean <- function(data, response,
                                        den_var = "Den",
                                        group_var = "Tactic",
                                        group_levels = c("Territorial", "Sneaker"),
                                        n_perm = 9999) {
  
  g1 <- group_levels[1]
  g2 <- group_levels[2]
  
  
  dat <- data
  
  obs_stat <- den_delta_mean_stat(dat, response = response,
                                  den_var = den_var, group_var = group_var,
                                  g1 = g1, g2 = g2)
  
  
  dens <- unique(dat[[den_var]])
  den_indices <- lapply(dens, function(d) which(dat[[den_var]] == d))
  names(den_indices) <- dens
  
  perm_stats <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    dat_b <- dat
    
    t_perm <- dat_b[[group_var]]
    
    
    for (d in dens) {
      idx <- den_indices[[d]]
      
      if(length(idx) > 1) {
        t_perm[idx] <- sample(t_perm[idx], length(idx), replace = FALSE)
      }
    }
    
    dat_b[[group_var]] <- t_perm
    
    perm_stats[b] <- den_delta_mean_stat(dat_b, response = response,
                                         den_var = den_var, group_var = group_var,
                                         g1 = g1, g2 = g2)
  }
  
  
  perm_stats <- perm_stats[!is.na(perm_stats)]
  
  
  if (is.na(obs_stat)) {
    p_two_sided <- NA
  } else {
    
    p_two_sided <- (sum(abs(perm_stats) >= abs(obs_stat)) + 1) / (length(perm_stats) + 1)
  }
  
  # --- 绘制置换分布直方图 ---
  # 创建数据框用于 ggplot
  plot_df <- data.frame(perm_stats = perm_stats)
  
  p <- ggplot(plot_df, aes(x = perm_stats)) +
    geom_histogram(binwidth = diff(range(perm_stats)) / 30, fill = "lightblue", color = "black", alpha = 0.7) + # 直方图
    geom_vline(xintercept = obs_stat, color = "red", linetype = "dashed", linewidth = 1) + # 观测统计量 (使用 linewidth 替代 size)
    annotate("text", x = obs_stat, y = Inf, label = paste0("Observed = ", round(obs_stat, 3)), # 使用 annotate 替代 geom_text
             color = "red", vjust = 1.5, hjust = ifelse(obs_stat > mean(perm_stats), 1.1, -0.1)) + # 标注观测值
    labs(title = paste("Permutation Distribution for", response), # 图表标题
         x = paste("Permuted Delta Mean (", g1, " - ", g2, ")"),  # X轴标签
         y = "Frequency") + # Y轴标签
    theme_minimal() + # 简洁主题
    theme(plot.title = element_text(hjust = 0.5)) # 标题居中
  
  print(p) # 打印图表，在RStudio或其他绘图设备中显示
  
  # --- 结束绘制 ---
  
  list(
    response = response,
    obs_stat = obs_stat,
    p_value = p_two_sided,
    perm_stats = perm_stats # 返回置换统计量，以防后续需要
  )
}

## ==========对三个响应变量分别检验 ========== 
results <- lapply(resp_list, function(r) {
  perm_test_by_den_delta_mean(df2, response = r, n_perm = n_perm)
})
write.csv(
  results,
  file.path( "置换检验_results.csv"), 
  row.names = FALSE,                             
  fileEncoding = "UTF-8"                          
)
res_df <- data.frame(
  response = sapply(results, `[[`, "response"),
  obs_stat = sapply(results, `[[`, "obs_stat"),
  p_value  = sapply(results, `[[`, "p_value"),
  row.names = NULL
)

# Holm 校正（针对3个指标）
res_df$p_value_holm <- p.adjust(res_df$p_value, method = "holm")

print(res_df)
write.csv(
  res_df,
  file.path( "置换检验_results2.csv"), # 保存模型结果到 CSV 文件
  row.names = FALSE,                              # 不写入行名
  fileEncoding = "UTF-8"                          # 指定文件编码
)
## ==========更直观输出 ========== 
for (i in seq_len(nrow(res_df))) {
  cat(
    "
", res_df$response[i], "
",
"obs Δ_mean (Den-level mean of (Territorial - Sneaker)) = ", round(res_df$obs_stat[i], 4), "
",
"permutation p (two-sided) = ", signif(res_df$p_value[i], 4), "
",
"Holm-adjusted p = ", signif(res_df$p_value_holm[i], 4), "
"
  )
}

