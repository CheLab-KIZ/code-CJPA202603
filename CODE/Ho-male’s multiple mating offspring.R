library(glmmTMB)
data_text <- "
famale_ID	Dens	Mono_mult	Breeding_season	Stream	Sampled	Mating_males	Ho	He	F
M5	D10	mono	1	1	8	1	0.844	0.645	-0.312
M9	D10	mono	1	1	5	1	0.833	0.618	-0.292
M11	D12	mono	1	1	24	1	0.865	0.672	-0.286
M2	D8	mult	1	1	61	2	0.792	0.698	-0.131
M3	D9	mult	1	1	37	2	0.835	0.663	-0.233
M4	D10	mult	1	1	37	3	0.898	0.759	-0.178
M6	D10	mult	1	1	49	2	0.906	0.729	-0.245
M7	D10	mult	1	1	34	3	0.838	0.722	-0.169
M10	D10	mult	1	1	76	2	0.842	0.708	-0.193
M12	D13	mono	1	2	9	1	0.667	0.598	-0.084
M1	D14	mono	2	1	4	1	0.896	0.578	-0.553
M14	D14	mono	2	1	3	1	0.833	0.556	-0.508
M5	D15	mono	2	1	4	1	0.875	0.667	-0.318
M6	D15	mult	2	1	40	3	0.872	0.735	-0.190
M10	D15	mult	2	1	51	2	0.804	0.656	-0.226
M11	D15	mult	2	1	22	3	0.773	0.659	-0.157
M13	D17	mono	2	2	16	1	0.788	0.606	-0.295
M19	D17	mono	2	2	6	1	0.875	0.638	-0.393
M15	D16	mult	2	2	54	3	0.895	0.721	-0.241
M16	D16	mult	2	2	4	3	0.854	0.638	-0.347
M17	D16	mult	2	2	17	3	0.711	0.653	-0.088"
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
m_full <- glmmTMB(Ho ~ Mono_mult + Sampled + (1 | Dens) + (1 | famale_ID), 
                  data = data, 
                  family = beta_family(link = "logit"))

# 模型 2：仅包含洞穴随机效应
m_dens <- glmmTMB(Ho ~ Mono_mult + Sampled + (1 | Dens), 
                  data = data, 
                  family = beta_family(link = "logit"))

# 模型 3：仅包含雌性编号随机效应
m_female <- glmmTMB(Ho ~ Mono_mult + Sampled + (1 | famale_ID), 
                    data = data, 
                    family = beta_family(link = "logit"))

# 模型 4：无随机效应的 Beta-GLM 模型
m_glm <- glmmTMB(Ho ~ Mono_mult + Sampled, 
                 data = data, 
                 family = beta_family(link = "logit"))

# 模型 5：无随机效应的 Beta-GLM 模型
m_glm_nosampled <- glmmTMB(Ho ~ Mono_mult, 
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

# -------------------------------------------------------------
# 5.拟合最优模型，并进行模型诊断
# -------------------------------------------------------------
library(glmmTMB)
library(DHARMa)
library(performance)

# 1. 拟合最优模型
m_glm_nosampled <- glmmTMB(Ho ~ Mono_mult, 
                           data = data, 
                           family = beta_family(link = "logit"))

cat("\n=================== 似然比检验 (LRT) ===================\n")
print(drop1(m_glm_nosampled, test = "Chisq"))

# 2. 评估拟合优度 (Pseudo-R2)
print(r2(m_glm_nosampled))

# 3. DHARMa 模拟残差诊断 (注意：参数名是 n，不是 nsim)
sim_res <- simulateResiduals(fittedModel = m_glm_nosampled, n = 1000)

# 绘制诊断图 (左侧QQ图看残差分布，右侧残差对预测值看异方差)
plot(sim_res)

# 4. 正式统计检验
cat("\n--- 残差均匀性检验 ---\n")
print(testUniformity(sim_res))

cat("\n--- 过度离散检验 ---\n")
print(testDispersion(sim_res))

cat("\n--- 模拟异常值检验 (Outlier Test) ---\n")
print(testOutliers(sim_res))
