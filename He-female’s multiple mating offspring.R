
library(glmmTMB)

# 1. 导入数据
data_text <- "
famale_ID	Dens	Mono_mult	Breeding_season	Stream	Sampled	Mating_males	Ho	He	F
F2	D8	mono	1	1	28	1	0.690	0.591	-0.127
F3	D8	mono	1	1	33	1	0.879	0.672	-0.301
F9	D11	mono	1	1	35	1	0.862	0.683	-0.261
F4	D9	mono	1	1	36	1	0.830	0.656	-0.237
F10	D12	mono	1	1	24	1	0.865	0.672	-0.286
F15	D16	mono	2	2	47	1	0.904	0.705	-0.282
F5	D10	mult	1	1	74	5	0.883	0.774	-0.145
F6	D10	mult	1	1	76	4	0.851	0.714	-0.190
F7	D10	mult	1	1	25	4	0.873	0.724	-0.176
F11	D13	mult	1	2	10	2	0.675	0.620	-0.067
F12	D14	mult	2	1	7	2	0.869	0.595	-0.445
F14	D15	mult	2	1	13	2	0.897	0.688	-0.305
F6	D15	mult	2	1	64	3	0.842	0.693	-0.214
F2	D15	mult	2	1	40	4	0.771	0.723	-0.058
F16	D16	mult	2	2	8	2	0.854	0.670	-0.268
F17	D16	mult	2	2	13	3	0.750	0.642	-0.183
F18	D16	mult	2	2	7	3	0.679	0.637	-0.047
F19	D17	mult	2	2	23	3	0.817	0.686	-0.195"

data <- read.table(text = data_text, header = TRUE, stringsAsFactors = FALSE)

data$famale_ID       <- as.factor(data$famale_ID)
data$Dens            <- as.factor(data$Dens)
data$Mono_mult       <- as.factor(data$Mono_mult)
data$Breeding_season <- as.factor(data$Breeding_season)
data$Stream          <- as.factor(data$Stream)

# -------------------------------------------------------------
# 2. 构建候选模型 
# -------------------------------------------------------------

# 模型 1：全随机效应模型 (同时包含洞穴 Dens 与 雌性编号 famale_ID)
m_full <- glmmTMB(He ~ Mono_mult + Sampled + (1 | Dens) + (1 | famale_ID), 
                  data = data, 
                  family = beta_family(link = "logit"))

# 模型 2：仅包含洞穴随机效应
m_dens <- glmmTMB(He ~ Mono_mult + Sampled + (1 | Dens), 
                  data = data, 
                  family = beta_family(link = "logit"))

# 模型 3：仅包含雌性编号随机效应
m_female <- glmmTMB(He ~ Mono_mult + Sampled + (1 | famale_ID), 
                    data = data, 
                    family = beta_family(link = "logit"))

# 模型 4：无随机效应的 Beta-GLM 模型
m_glm <- glmmTMB(He ~ Mono_mult + Sampled, 
                 data = data, 
                 family = beta_family(link = "logit"))

# 模型 5：无随机效应的 Beta-GLM 模型
m_glm_nosampled <- glmmTMB(He ~ Mono_mult, 
                 data = data, 
                 family = beta_family(link = "logit"))
# -------------------------------------------------------------
# 3. 自定义小样本 AICc 计算与模型比较
# -------------------------------------------------------------
calculate_aicc <- function(model) {
  aic_val <- AIC(model)
  k <- as.numeric(attr(logLik(model), "df"))
  n <- nobs(model)
  aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
  return(aicc_val)
}

models <- list(
  "GLMM (Dens + Female)" = m_full,
  "GLMM (仅 Dens)"       = m_dens,
  "GLMM (仅 Female)"     = m_female,
  "Beta-GLM (无随机效应)" = m_glm,
  "Beta-GLM (无样本量考虑)" = m_glm_nosampled
)

aicc_values <- sapply(models, calculate_aicc)
delta_aicc <- aicc_values - min(aicc_values)

aicc_table <- data.frame(
  Model = names(models),
  K = sapply(models, function(m) as.numeric(attr(logLik(m), "df"))),
  AICc = round(aicc_values, 2),
  Delta_AICc = round(delta_aicc, 2)
)

# 按 Delta AICc 升序排序
aicc_table <- aicc_table[order(aicc_table$AICc), ]

cat("\n===候选模型的 AICc 比较结果 ===\n")
print(aicc_table)

# -------------------------------------------------------------
# 4. 输出 AICc 最优模型的详细摘要
# -------------------------------------------------------------
best_model_name <- rownames(aicc_table)[1]
cat(paste("\n=== 最优模型详细结果：", best_model_name, "===\n"))
print(summary(models[[best_model_name]]))

#检验剔除各个固定效应后模型的卡方值和显著性变化
drop1(m_glm, test = "Chisq")

# -------------------------------------------------------------
# 5.拟合最优模型，并进行模型诊断
# -------------------------------------------------------------

library(glmmTMB)
library(DHARMa)
library(performance)

# 1. 拟合最优模型
m_optimal <- glmmTMB(He ~ Mono_mult + Sampled, 
                     data = data, 
                     family = beta_family(link = "logit"))

# 2. 评估拟合优度 (Pseudo-R2)
print(r2(m_optimal))

# 3. DHARMa 模拟残差诊断 (注意：参数名是 n，不是 nsim)
sim_res <- simulateResiduals(fittedModel = m_optimal, n = 1000)

# 绘制诊断图 (左侧QQ图看残差分布，右侧残差对预测值看异方差)
plot(sim_res)

# 4. 正式统计检验
cat("\n--- 残差均匀性检验 ---\n")
print(testUniformity(sim_res))

cat("\n--- 过度离散检验 ---\n")
print(testDispersion(sim_res))

cat("\n--- 模拟异常值检验 (Outlier Test) ---\n")
print(testOutliers(sim_res))


# -------------------------------------------------------------
# 6.检验上面2个固定效应分别的解释率
# -------------------------------------------------------------

library(glmmTMB)
library(performance)

# 1. 拟合三个模型
m_full      <- glmmTMB(He ~ Mono_mult + Sampled, data = data, family = beta_family(link = "logit"))
m_no_mono   <- glmmTMB(He ~ Sampled, data = data, family = beta_family(link = "logit"))
m_no_sampled<- glmmTMB(He ~ Mono_mult, data = data, family = beta_family(link = "logit"))

# 2. 稳妥提取 R2 数值（强制转换为数值型）
r2_full       <- as.numeric(r2(m_full))
r2_no_mono    <- as.numeric(r2(m_no_mono))
r2_no_sampled <- as.numeric(r2(m_no_sampled))

# 3. 计算独立解释力与重叠部分
unique_mono    <- r2_full - r2_no_mono
unique_sampled <- r2_full - r2_no_sampled
shared_r2      <- r2_full - unique_mono - unique_sampled

# 4. 打印结果
cat(sprintf("全模型总 R2: %.3f\n", r2_full))

cat(sprintf("交配系统 (Mono_mult) 独占解释的 R2: %.3f (约占总 R2 的 %.1f%%)\n", unique_mono, (unique_mono/r2_full)*100))

cat(sprintf("子代样本数 (Sampled) 独占解释的 R2: %.3f (约占总 R2 的 %.1f%%)\n", unique_sampled, (unique_sampled/r2_full)*100))

cat(sprintf("两者的重叠/协同解释 R2: %.3f\n", shared_r2))


