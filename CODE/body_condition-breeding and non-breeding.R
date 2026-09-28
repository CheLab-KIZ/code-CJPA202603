
# ==============================================================================
# 步骤一：响应变量尺度与正态性评估脚本
# ==============================================================================

# 1. 加载必要的包
library(ggplot2)
library(MASS)      # 提供 boxcox() 函数，用于寻找最佳变换指数

# 2. 读取并清洗数据
file_path <- "row_data-breeding_and_non-breeding.csv"
data <- read.csv(file_path, header = TRUE, fileEncoding = "UTF-8-BOM")

data$Breeding_tactics <- as.factor(data$Breeding_tactics)
data$Body_condition      <- as.numeric(data$Body_condition)

# 剔除 NA
data_clean <- na.omit(data[, c("Body_condition", "Breeding_tactics")])
cat("有效样本量 (N):", nrow(data_clean), "\n\n")

#【方法一】Box-Cox 变换分析（核心：寻找最佳 Lambda）
# 拟合一个基础线性模型来执行 Box-Cox
base_lm <- lm(Body_condition ~ Breeding_tactics, data = data_clean)

par(mfrow = c(1, 2))
# 绘制 Box-Cox 曲线
boxcox_res <- boxcox(base_lm, plotit = TRUE, 
                     main = "Box-Cox Transformation Curve")

# 提取使对数似然值达到最大时的 lambda 值
lambda <- boxcox_res$x[which.max(boxcox_res$y)]
cat("==================================================\n")
cat("【Box-Cox 建议】最优 Lambda (λ) 估计值为:", round(lambda, 3), "\n")
if (abs(lambda) < 0.2) {
  cat("结论建议: Lambda 接近 0，强烈建议进行【对数转换 (log)】。\n")
} else if (abs(lambda - 1) < 0.2) {
  cat("结论建议: Lambda 接近 1，数据无需转换，可直接使用原始数据。\n")
} else {
  cat("结论建议: Lambda 接近", round(lambda, 2), "，可考虑相应的幂转换。\n")
}
cat("==================================================\n\n")

#【方法二】正态性检验对比 (Shapiro-Wilk Test)
cat("--- 原始数据 Shapiro-Wilk 正态性检验 ---\n")
print(shapiro.test(data_clean$Body_condition))

data_clean$Log_Total_length <- log(data_clean$Body_condition)

cat("\n--- 对数转换后 Shapiro-Wilk 正态性检验 ---\n")
print(shapiro.test(data_clean$Log_Body_condition))


# ==============================================================================
# 步骤二：AICc 模型选择
# ==============================================================================

library(lme4)
library(lmerTest)

# 1. 读取数据并确保分类变量是因子
file_path <- "row_data-breeding_and_non-breeding.csv"
data <- read.csv(file_path, header = TRUE, fileEncoding = "UTF-8-BOM")

data$Breeding_tactics <- as.factor(data$Breeding_tactics)
data$Stream           <- as.factor(data$Stream)
data$Season           <- as.factor(data$Season) # 补充 Season 因子化
data$Body_condition     <- as.numeric(data$Body_condition)

# 2. 剔除 NA 并构建干净的数据框 (必须先清洗，再基于清洗后的数据做对数转换)
data_clean <- na.omit(data[, c("Body_condition", "Breeding_tactics", "Stream", "Season")])



# 3. 自定义 AICc 计算函数
calculate_AICc <- function(model) {
  k <- attr(logLik(model), "df")
  n <- nobs(model)
  aic_val <- AIC(model)
  aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
  return(aicc_val)
}

# 4. 构建候选模型（统一使用 data_clean，变量使用 Body_condition）
# 模型1：仅 Stream 随机效应
m_stream <- lmer(Body_condition ~ Breeding_tactics + (1 | Stream), data = data_clean, REML = TRUE)


# 模型3：无随机效应的普通线性模型 (LM)
m_lm <- lm(Body_condition ~ Breeding_tactics, data = data_clean)

# 5. 提取各模型的 AICc 并汇总对比
results <- data.frame(
  Model = c("Stream only",  "No random effect (LM)"),
  AICc = c(calculate_AICc(m_stream),calculate_AICc(m_lm))
)

# 按 AICc 从小到大排序（AICc 越小模型越优）
results <- results[order(results$AICc), ]

print(results)

# ==============================================================================
# 步骤三：最优模型模型诊断
# ==============================================================================
# 1. 加载必要的包
library(lme4)
library(lmerTest)  
library(ggplot2)

# 3. 拟合对数转换后的最优 GLMM 模型 (使用 REML)
optimal_model <- lm(Body_condition ~ Breeding_tactics, 
                      data = data_clean,
)

# 4. 输出模型统计结果摘要
cat("\n=================== 模型统计摘要 (Summary) ===================\n")
print(summary(optimal_model))

# 5. 模型方差分析表（Type III ANOVA 检验固定效应的显著性）
cat("\n=================== 方差分析表 (ANOVA) ===================\n")
print(anova(optimal_model))

# 6. 计算模型解释力 (R²)，替代 MuMIn
if (requireNamespace("performance", quietly = TRUE)) {
  library(performance)
  cat("\n=================== 模型 R² 评估 (performance包) ===================\n")
  print(r2(optimal_model))
} else {
  cat("\n提示: 未安装 performance 包，如需 R² 可通过手动计算或报告随机效应标准差。\n")
}

# 7. 模型残差诊断（基于对数转换后的模型）
par(mfrow = c(1, 2))
plot(fitted(optimal_model), residuals(optimal_model), 
     xlab = "Fitted Values (Log-scale)", ylab = "Residuals", 
     main = "Residuals vs Fitted (Log-transformed)")
abline(h = 0, col = "red", lty = 2)

# Q-Q 图（检查残差正态性）
qqnorm(residuals(optimal_model))
qqline(residuals(optimal_model), col = "red")
par(mfrow = c(1, 1))


# ==============================================================================
# 步骤四：绘制原始数据箱线图
# ==============================================================================

cat("
=================== 原始数据箱线图 (Body_condition) ===================
")


custom_colors <- c("Breeding" = "#F8B6B2", "Non-breeding" = "#B2D9D9") 


p <- ggplot(data_clean, aes(x = Breeding_tactics, y = Body_condition, fill = Breeding_tactics)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.8, color = "black", size = 2) +
  scale_fill_manual(values = custom_colors) + 
  labs(x = "Breeding Tactics", y = "Body_condition") + 
  ggtitle("Body_condition by Breeding Tactics") + 
  theme_minimal() + 
  theme(
    legend.position = "none", 
    axis.line = element_line(colour = "black") 
  )
ggsave("Body_condition_boxplot.svg", plot = p, width = 6, height = 4, dpi = 300)

