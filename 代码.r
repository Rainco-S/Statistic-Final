library(MASS)
library(lmtest)
library(lubridate)
library(ggplot2)
library(gganimate)
library(forecast)
library(tseries)
library(tidyverse)
library(imputeTS)
library(randomForest)
library(zoo)
library(dplyr)
library(purrr)
library(lubridate)
library(car)
library(forecast)
library(corrplot)

path <- '/Users/DELL/Desktop/普通统计学/期末作业/'
data <- read.csv(file.path(path,"/Beijing_Wanliu_data.csv"), stringsAsFactors = FALSE)

str(data)

# 计算缺失比例
missing_ratio <- sapply(data, function(x) mean(is.na(x)) * 100)
missing_df <- data.frame(缺失比例 = round(missing_ratio, 3))
print(missing_df)

# 缺失值处理
numeric_vars <- c("PM2.5", "PM10", "SO2", "CO", "NO2", "O3", "TEMP", "DEWP", "HUMI", "PRES", "WSPM", "year", "month", "day", "hour")
categorical_vars <- c("wd")
# 数值型变量：线性插值
data_numeric <- data[, numeric_vars]
data_numeric_imputed <- as.data.frame(
  apply(data_numeric, 2, function(col) {
    na.interp(col)
  })
)

# 分类型变量：最近邻填充
data_categorical <- data[, categorical_vars, drop = FALSE]
data_categorical_imputed <- na_locf(data_categorical)  # 向前填充
data_categorical_imputed <- na_locf(data_categorical_imputed, option = "locf_back")  # 向后填充

# 合并处理后的数据
data_imputed <- cbind(data_numeric_imputed, data_categorical_imputed)
data_imputed <- data_imputed %>%
  mutate(
    season = case_when(month %in% 3:5 ~ "春季",
                       month %in% 6:8 ~ "夏季",
                       month %in% 9:11 ~ "秋季",
                       TRUE ~ "冬季"),
    period = ifelse(hour %in% 6:18, "day", "night"),
    pm25level = case_when(PM2.5 <= 35 ~ "优良",
                           PM2.5 <= 75 ~ "轻度污染",
                           PM2.5 <= 115 ~ "中度污染",
                           TRUE ~ "重度污染"),
    date = as.Date(paste(year, month, day, sep = "-")),
  )
pm25_daily <- data_imputed %>% group_by(date) %>% summarise(PM2.5 = mean(PM2.5))

# PM2.5浓度分布直方图/频数分布直方图
hist(data_imputed$PM2.5,
     main = "PM2.5浓度频数分布直方图",
     xlab = "PM2.5浓度（μg/m³）",
     ylab = "频数",
     col = "lightblue",
     breaks = 20)
hist_info <- hist(data_imputed$PM2.5, plot = FALSE, breaks = 20)
text(hist_info$mids, hist_info$counts, labels = hist_info$counts, 
     pos = 3, cex = 0.7)

# 不同季节PM2.5均值对比条形图
season_pm25 <- data_imputed %>%
  group_by(season) %>%
  summarise(pm25_mean = mean(PM2.5), .groups = "drop")
barplot(season_pm25$pm25_mean,
        names.arg = season_pm25$season,
        main = "不同季节PM2.5均值条形图",
        xlab = "季节",
        ylab = "PM2.5均值（μg/m³）",
        col = c("green", "red", "orange", "blue"))
text(1:4, season_pm25$pm25_mean + 5, 
     labels = round(season_pm25$pm25_mean, 1))

# PM2.5污染等级占比饼图
pm25_level_count <- table(data_imputed$pm25level)
pie(pm25_level_count,
    labels = paste(names(pm25_level_count), "\n", round(pm25_level_count/sum(pm25_level_count)*100, 1), "%"),
    main = "PM2.5污染等级占比饼图",
    col = c("lightgreen", "yellow", "orange", "red"))

# 季节×污染等级交叉分布列联表
season_level_table <- table(data_imputed$season, data_imputed$pm25level)
print(season_level_table)
mosaicplot(season_level_table,
           main = "季节×PM2.5污染等级列联表（马赛克图）",
           xlab = "季节",
           ylab = "污染等级",
           col = c("lightgreen", "yellow", "orange", "red"))

# PM2.5与气温的相关性散点图
plot(data_imputed$TEMP, data_imputed$PM2.5,
     main = "PM2.5浓度与气温散点图",
     xlab = "气温（℃）",
     ylab = "PM2.5浓度（μg/m³）",
     col = alpha("darkblue", 0.3),
     pch = 16)
abline(lm(PM2.5 ~ TEMP, data_imputed), col = "red", lwd = 2)

# PM2.5日均值时间趋势时间序列图
pm25_daily <- data_imputed %>%
  group_by(date) %>%
  summarise(pm25_daily = mean(PM2.5), .groups = "drop")
plot(pm25_daily$date, pm25_daily$pm25_daily,
     type = "l",
     main = "PM2.5日均值时间序列图",
     xlab = "日期",
     ylab = "PM2.5日均值（μg/m³）",
     col = "darkred",
     lwd = 1)
lines(loess.smooth(pm25_daily$date, pm25_daily$pm25_daily), col = "blue", lwd = 2)

# 风向主次图（PM2.5均值主图+频次次轴）
wd_stats <- data_imputed %>%
  group_by(wd) %>%
  summarise(pm25_mean = mean(PM2.5),
            wd_count = n(),
            .groups = "drop") %>%
  arrange(desc(pm25_mean)) %>%
  head(8)

par(mar = c(5, 4, 4, 4) + 0.1)
barplot(wd_stats$pm25_mean,
        names.arg = wd_stats$wd,
        main = "风向-PM2.5均值主次图",
        xlab = "风向",
        ylab = "PM2.5均值（μg/m³）",
        col = "lightcoral")
par(new = TRUE)
plot(1:8, wd_stats$wd_count,
     type = "l",
     col = "darkblue",
     lwd = 2,
     axes = FALSE,
     xlab = "", ylab = "")
axis(4, family = "PingFang")
mtext("风向频次（次）", side = 4, line = 3)
legend("topright",
       legend = c("PM2.5均值（主轴）", "风向频次（次轴）"),
       col = c("lightcoral", "darkblue"),
       lty = c(NA, 1), pch = c(15, NA))

# 数值型变量描述性统计
numeric_stats <- data_imputed %>% 
  summarise(across(all_of(numeric_vars),
                   list(均值 = mean,
                        中位数 = median,
                        方差 = var,
                        标准差 = sd),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(everything(),
               names_to = c("变量", "统计量"),
               names_sep = "_") %>%
  pivot_wider(names_from = 统计量, values_from = value) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2)))

print(numeric_stats)

# 数值型变量箱线图
par(mfrow = c(2, 5))
for (var in c("PM2.5", "PM10", "SO2", "CO", "NO2", "O3", "TEMP", "HUMI", "PRES", "WSPM")) {
  boxplot(data_imputed[[var]], main = var, col = "lightgreen")
}

# 分类型变量汇总
categorical_vars_new <- c("wd", "season", "period", "pm25level")

categorical_stats <- map_dfr(categorical_vars_new, function(var) {
  freq <- as.vector(table(data_imputed[[var]]))          # 转成纯向量
  cat_names <- names(table(data_imputed[[var]]))        # 类别名
  prop <- round(prop.table(table(data_imputed[[var]])) * 100, 2)
  mode_val <- cat_names[which.max(freq)]
  
  tibble(
    变量 = var,
    类别 = cat_names,
    频数 = freq,
    比例 = prop,
    众数 = mode_val
  )
})

print(categorical_stats)

# 验证新变量
table(data_imputed$season)
table(data_imputed$period)

# 假设检验与统计推断
# PM2.5 Q-Q图
qqnorm(data_imputed$PM2.5, main = "PM2.5 Q-Q图")
qqline(data_imputed$PM2.5, col = "red", lwd = 2)
# TEMP Q-Q图
qqnorm(data_imputed$TEMP, main = "气温 Q-Q图")
qqline(data_imputed$TEMP, col = "red", lwd = 2)
# WSPM Q-Q图
qqnorm(data_imputed$WSPM, main = "风速 Q-Q图")
qqline(data_imputed$WSPM, col = "red", lwd = 2)

# 计算样本统计量
pm25_mean <- mean(data_imputed$PM2.5, na.rm = TRUE)
pm25_sd <- sd(data_imputed$PM2.5, na.rm = TRUE)
n <- nrow(data_imputed)
z_alpha2 <- qnorm(0.975)  # 95%置信水平

# 置信区间
margin_error <- z_alpha2 * pm25_sd / sqrt(n)
ci_known_var <- c(pm25_mean - margin_error, pm25_mean + margin_error)
print(round(ci_known_var, 2))

# 筛选冬季数据
pm25_winter <- data_imputed %>% filter(season == "冬季") %>% pull(PM2.5)

# t置信区间
winter_mean <- mean(pm25_winter)
winter_sd <- sd(pm25_winter)
t_ci <- t.test(pm25_winter)$conf.int
print(round(t_ci, 2))

# 对比全年与冬季均值
cat("全年PM2.5均值：", round(pm25_mean, 2), "\n")
cat("冬季PM2.5均值：", round(mean(pm25_winter), 2), "\n")

# 比例区间估计
p1_hat <- mean(data_imputed$PM2.5 > 75, na.rm = TRUE)
p2_hat <- mean(data_imputed$PM2.5 <= 35, na.rm = TRUE)

# 置信区间
ci_p1 <- p1_hat + c(-1, 1) * z_alpha2 * sqrt(p1_hat*(1-p1_hat)/n)
ci_p2 <- p2_hat + c(-1, 1) * z_alpha2 * sqrt(p2_hat*(1-p2_hat)/n)

cat("重污染比例p₁：", round(p1_hat*100, 1), "%, 95%置信区间：", round(ci_p1*100, 1), "\n")
cat("优良空气比例p₂：", round(p2_hat*100, 1), "%, 95%置信区间：", round(ci_p2*100, 1), "\n")

# 样本量估计
E <- 0.05  # 误差控制在5%
n1 <- (z_alpha2^2 * p1_hat * (1-p1_hat)) / E^2
n2 <- (z_alpha2^2 * p2_hat * (1-p2_hat)) / E^2
n_conservative <- (z_alpha2^2 * 0.5 * 0.5) / E^2

cat("重污染比例所需样本量n1：", ceiling(n1), "\n")
cat("优良空气比例所需样本量n2：", ceiling(n2), "\n")
cat("最保守情形所需样本量n：", ceiling(n_conservative), "\n")

# 分组数据
pm25_day <- data_imputed %>% filter(period == "day") %>% pull(PM2.5)
pm25_night <- data_imputed %>% filter(period == "night") %>% pull(PM2.5)

# F检验方差齐性
pm25_day_var <- var(pm25_day)
pm25_night_var <- var(pm25_night)
var_test <- var.test(pm25_night, pm25_day)
cat("F检验结果：F=", round(var_test$statistic, 3), ", p值=", var_test$p.value, "\n")

# t检验（根据F检验结果选方差不齐）
t_test <- t.test(pm25_day, pm25_night, var.equal = FALSE)
cat("t检验结果：t=", round(t_test$statistic, 3), ", p值=", t_test$p.value, "\n")
cat("白天PM2.5均值：", round(mean(pm25_day), 1), "μg/m³\n")
cat("夜间PM2.5均值：", round(mean(pm25_night), 1), "μg/m³\n")

# 单因素ANOVA
anova_model <- aov(PM2.5 ~ season, data = data_imputed)
anova_table <- summary(anova_model)[[1]]

# 计算总变异 (Total SS = Between SS + Residual SS)
total_ss <- sum(anova_table$`Sum Sq`)
total_df <- sum(anova_table$Df)

f_value <- anova_table$`F value`[1]
p_value <- anova_table$`Pr(>F)`[1]

anova_table_clean <- data.frame(
  来源 = c("组间", "组内", "总和"),
  自由度 = c(anova_table$Df, total_df),
  平方和 = c(
    round(anova_table$`Sum Sq`, 3),
    round(total_ss, 3)
  ),
  均方 = c(
    round(anova_table$`Mean Sq`, 3),
    ""    # ← 总变异这一行不需要均方
  ),
  F值 = c(round(f_value, 3), "", ""),
  p值 = c(
    ifelse(p_value < 0.001, "<0.001", round(p_value, 4)),
    "", ""
  )
)
print(anova_table_clean)

# 定义自变量（数值型变量排除时间变量）
predictors <- c("PM10", "SO2", "CO", "NO2", "O3", "TEMP", "DEWP", "HUMI", "PRES", "WSPM")
lm_model <- lm(PM2.5 ~ ., data = data_imputed[, c("PM2.5", predictors)])

# 回归输出
summary_lm <- summary(lm_model)
print(summary_lm)
conf_int <- confint(lm_model, level = 0.95)
print(conf_int)

# 残差诊断图
par(mfrow = c(2, 2))
plot(lm_model)

# 逐步回归（双向）
lm_null <- lm(PM2.5 ~ 1, data = data_imputed)
step_model <- stepAIC(lm_model,
                  scope = list(lower = lm_null, upper = lm_model),
                  direction = "both",
                  trace = TRUE,
                  test = "F")   # 基于 F 值

# 提取保留变量
summary(step_model)
summary_step <- summary(step_model)
retained_vars <- names(coef(step_model))[-1]  # 排除截距项
cat("逐步回归保留变量：", paste(retained_vars, collapse = ", "), "\n")

cat("全模型 R²：", round(summary_lm$r.squared, 3), ", 调整后R²：", round(summary_lm$adj.r.squared, 3), ", AIC：", round(AIC(lm_model), 0), ", BIC：", round(BIC(lm_model), 0), "\n")
cat("逐步模型 R²：", round(summary_step$r.squared, 3), ", 调整后R²：", round(summary_step$adj.r.squared, 3), ", AIC：", round(AIC(step_model), 0), ", BIC：", round(BIC(lm_model), 0), "\n")

# 二次项和滞后项的回归模型
# 生成二阶项
data_ext <- data_imputed %>%
  mutate(
    TEMP2 = TEMP^2,
    HUMI2 = HUMI^2,
    WSPM2 = WSPM^2,
    PRES2 = PRES^2,
    DEWP2 = DEWP^2
  )

# 生成一阶滞后项
vars_to_lag <- c("PM2.5", "TEMP", "DEWP", "HUMI", "PRES", "WSPM", "PM10", "SO2", "CO", "NO2", "O3")

for (v in vars_to_lag) {
  lag_name <- paste0(v, "_lag1")
  data_ext[[lag_name]] <- dplyr::lag(data_ext[[v]], 1)
}

# 去掉第一行（滞后项为空）
data_ext <- data_ext %>% filter(!is.na(PM2.5_lag1))

quad_predictors <- c("TEMP2", "DEWP2", "HUMI2", "PRES2", "WSPM2")
lag_predictors <- paste0(vars_to_lag, "_lag1")
all_predictors <- c(predictors, quad_predictors, lag_predictors)

lm_ext <- lm(PM2.5 ~ ., data = data_ext[, c("PM2.5", all_predictors)])
summary(lm_ext)
summary_ext <- summary(lm_ext)

lm_ext_null <- lm(PM2.5 ~ 1, data = data_ext[, c("PM2.5", all_predictors)])
step_ext <- stepAIC(lm_ext,
                    scope = list(lower = lm_ext_null, upper = lm_ext),
                    direction = "both",
                    trace = TRUE,
                    test = "F")
summary_step_ext <- summary(step_ext)

cat("扩展逐步回归保留变量：", paste(names(coef(step_ext))[-1], collapse = ", "), "\n")

mse <- function(model, data) {
  mean(model$residuals^2)
}

mse_lm  <- mse(step_model, data_imputed)
mse_ext <- mse(lm_ext, data_ext)
mse_step_ext <- mse(step_ext, data_ext)

cat("线性模型 R²：", round(summary_lm$r.squared, 3),
    " 调整后 R²：", round(summary_lm$adj.r.squared, 3),
    " AIC：", round(AIC(step_model), 1),
    " BIC：", round(BIC(step_model), 1),
    " MSE：", round(mse_lm, 2), "\n")

cat("扩展模型 R²：", round(summary_ext$r.squared, 3),
    " 调整后 R²：", round(summary_ext$adj.r.squared, 3),
    " AIC：", round(AIC(lm_ext), 1),
    " BIC：", round(BIC(lm_ext), 1),
    " MSE：", round(mse_ext, 2), "\n")

cat("扩展逐步筛选后模型 R²：", round(summary_step_ext$r.squared, 3),
    " 调整后 R²：", round(summary_step_ext$adj.r.squared, 3),
    " AIC：", round(AIC(step_ext), 1),
    " BIC：", round(BIC(step_ext), 1),
    " MSE：", round(mse_step_ext, 2), "\n")

# 多阶滞后
data_ext_mul <- data_imputed %>%
  mutate(
    TEMP2 = TEMP^2,
    HUMI2 = HUMI^2,
    WSPM2 = WSPM^2,
    PRES2 = PRES^2,
    DEWP2 = DEWP^2
  ) %>%
  arrange(year, month, day, hour)

# 生成多阶滞后项
lag_orders <- c(1, 3, 6, 12)
for (v in vars_to_lag) {
  for (lag in lag_orders) {
    lag_name <- paste0(v, "_lag", lag)
    data_ext_mul[[lag_name]] <- dplyr::lag(data_ext_mul[[v]], lag)
  }
}
data_ext_mul <- data_ext_mul %>% filter(complete.cases(data_ext_mul[, grepl("lag", colnames(data_ext_mul))]))

lag_predictors_mul <- colnames(data_ext_mul)[grepl("lag[13612]$", colnames(data_ext_mul))]
all_predictors_mul <- c(predictors, quad_predictors, lag_predictors_mul)
lm_ext_mul <- lm(PM2.5 ~ ., data = data_ext_mul[, c("PM2.5", all_predictors_mul)])
summary(lm_ext_mul)
summary_ext_mul <- summary(lm_ext_mul)

# 多阶滞后项模型逐步回归
lm_ext_null_mul <- lm(PM2.5 ~ 1, data = data_ext_mul[, c("PM2.5", all_predictors_mul)])
step_ext_mul <- stepAIC(lm_ext_mul,
                    scope = list(lower = lm_ext_null_mul, upper = lm_ext_mul),
                    direction = "both",
                    trace = TRUE,
                    test = "F")
summary(step_ext_mul)
summary_step_ext_mul <- summary(step_ext_mul)

retained_lag_vars_mul <- grep("lag", names(coef(step_ext_mul)), value = TRUE)
cat("逐步回归保留的多阶滞后项：", paste(retained_lag_vars_mul, collapse = ", "), "\n")

mse <- function(model, data) {
  mean(model$residuals^2)
}

mse_lm  <- mse(step_model, data_imputed)
mse_ext_mul <- mse(lm_ext_mul, data_ext_mul)
mse_step_ext_mul <- mse(step_ext_mul, data_ext_mul)

cat("线性模型 R²：", round(summary_lm$r.squared, 3),
    " 调整后 R²：", round(summary_lm$adj.r.squared, 3),
    " AIC：", round(AIC(step_model), 1),
    " BIC：", round(BIC(step_model), 1),
    " MSE：", round(mse_lm, 2), "\n")

cat("筛选后一阶滞后模型 R²：", round(summary_step_ext$r.squared, 3),
    " 调整后 R²：", round(summary_step_ext$adj.r.squared, 3),
    " AIC：", round(AIC(step_ext), 1),
    " BIC：", round(BIC(step_ext), 1),
    " MSE：", round(mse_step_ext, 2), "\n")

cat("筛选后多阶滞后项模型 R²：", round(summary_step_ext_mul$r.squared, 3),
    " 调整后 R²：", round(summary_step_ext_mul$adj.r.squared, 3),
    " AIC：", round(AIC(step_ext_mul), 1),
    " BIC：", round(BIC(step_ext_mul), 1),
    " MSE：", round(mse_step_ext_mul, 2), "\n")

# 随机森林模型
rf_model <- randomForest(PM2.5 ~ ., data = data_ext[, c("PM2.5", all_predictors)], ntree = 100, importance = TRUE)
rf_r2 <- 1 - rf_model$mse[rf_model$ntree]/var(data_ext$PM2.5, na.rm = TRUE)

# 变量重要性
var_importance <- importance(rf_model) %>% 
  as.data.frame() %>% 
  arrange(desc(IncNodePurity)) %>%
  head(10)

cat("随机森林模型调整后R²：", round(rf_r2, 2), "\n")
print(var_importance)

# 变量重要性图
varImpPlot(rf_model, main = "随机森林变量重要性（Top10）", n.var = 10)

# SARIMA/SARIMAX季节性时间序列模型
pm25_ts <- ts(pm25_daily$pm25_daily, frequency = 365)

# 平稳性检验ADF
adf_result <- adf.test(pm25_ts)
print(adf_result) # p<0.05说明数据平稳，无需差分

# 自动选阶+拟合SARIMA
sarima_model <- auto.arima(
  pm25_ts,
  seasonal = TRUE,
  stepwise = FALSE,
  approximation = FALSE
)

pred_sarima <- forecast(sarima_model, h = 30)

plot(pred_sarima,
     main = "PM2.5 日均值 SARIMA 模型预测",
     xlab = "日期",
     ylab = "PM2.5 (μg/m³)")

# 给pm25_daily加季节列
pm25_daily <- pm25_daily %>%
  mutate(date = ymd(date)) %>%
  mutate(season = case_when(
    month(date) %in% 3:5 ~ "春季",
    month(date) %in% 6:8 ~ "夏季",
    month(date) %in% 9:11 ~ "秋季",
    TRUE ~ "冬季"
  ))

p <- ggplot(pm25_daily, aes(x = date, y = pm25_daily, color = season)) +
  geom_line(lwd = 1) + # 线条加粗
  geom_point(size = 1.2, alpha = 0.8) + # 点标记
  scale_color_manual(values = c("春季"="green", "夏季"="blue", "秋季"="orange", "冬季"="red")) +
  labs(x = "日期", y = "PM2.5日均值（μg/m³）", title = "2014年PM2.5浓度动态变化") +
  theme_bw() +
  theme(text = element_text(family = "PingFang"), # macOS中文字体
        plot.title = element_text(hjust = 0.5)) +
  transition_reveal(date) # 核心：按日期滑动展示

# 导出GIF
animate(p, fps = 8, duration = 12, width = 900, height = 500)
anim_save(file.path(path,"figures/PM25动态变化.gif"))
