library(mice)
library(dplyr)
set.seed(123)
setwd("D:/")
df <- read.csv("D:/data.csv") 
data <- df


for (var in cont_vars) {
  x <- data[[var]]
  q25 <- quantile(x, 0.25, na.rm = TRUE)
  q75 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q75 - q25
  
  upper_limit <- q75 + 3 * iqr
  lower_limit <- q25 - 3 * iqr
  
  data[[var]][x > upper_limit | x < lower_limit] <- NA
}
na_proportion <- colMeans(is.na(data)) 
selected_vars <- names(na_proportion[na_proportion <= 0.3]) 
df <- data[, selected_vars]
cat_vars <- c("gender", "ethnicity", "age", "education", "smoking", "drinking")
cont_vars <- c("BMI", "HR", "SBP", "DBP", "ALP", "ALB", "BUN", "TC","HDL", "BG", "TG", "SCr", "Ccr", "WBC", "HBA1C")

df <- df %>%
  mutate(across(all_of(cat_vars), as.factor))
###########################################################
mice.impute.mode <- function(y, ry, x = NULL, ...) {
  y_obs <- y[ry]
  if (length(y_obs) == 0) return(rep(NA, sum(!ry)))
  tbl <- table(y_obs)
  mode_val <- names(tbl)[which.max(tbl)]
  if (is.factor(y)) mode_val <- factor(mode_val, levels = levels(y))
  rep(mode_val, sum(!ry))
}

meth <- rep("", ncol(df))
names(meth) <- colnames(df)
meth[cat_vars] <- "mode"
meth[cont_vars] <- "pmm"

tempData <- mice(
  data = df,
  m = 5,        
  maxit = 50,
  seed = 500,
  meth = meth) 
densityplot(tempData)  

data1 <- complete(tempData, 1) 
data2 <- complete(tempData, 2)  
data3 <- complete(tempData, 3)  
data4 <- complete(tempData, 4)  
data5 <- complete(tempData, 5)


fit1 <- glm(CAN~gender + ethnicity + age + education + smoking + drinking + BMI + HR + SBP + DBP + ALP + ALB + BUN + TC + HDL + BG + TG + SCr + Ccr + WBC + HBA1C , data = data1, family="binomial")
fit2 <- glm(CAN~gender + ethnicity + age + education + smoking + drinking + BMI + HR + SBP + DBP + ALP + ALB + BUN + TC + HDL + BG + TG + SCr + Ccr + WBC + HBA1C , data = data2, family="binomial") 
fit3 <- glm(CAN~gender + ethnicity + age + education + smoking + drinking + BMI + HR + SBP + DBP + ALP + ALB + BUN + TC + HDL + BG + TG + SCr + Ccr + WBC + HBA1C , data = data3, family="binomial")
fit4 <- glm(CAN~gender + ethnicity + age + education + smoking + drinking + BMI + HR + SBP + DBP + ALP + ALB + BUN + TC + HDL + BG + TG + SCr + Ccr + WBC + HBA1C , data = data4, family="binomial") 
fit5 <- glm(CAN~gender + ethnicity + age + education + smoking + drinking + BMI + HR + SBP + DBP + ALP + ALB + BUN + TC + HDL + BG + TG + SCr + Ccr + WBC + HBA1C , data = data5, family="binomial")

#pooled using the Rubin' Rule
glmresult <- summary(pool(list(fit1, fit2, fit3, fit4, fit5)), conf.int = T, conf.level = 0.95) 
orform <- list()
orform[["OR"]] <- exp(glmresult[["estimate"]])
orform[["LowerLimit"]] <- exp(glmresult[["conf.low"]])
orform[["UpperLimit"]] <- exp(glmresult[["conf.high"]])
orform


table(data1$CAN)
if(!require(DMwR)) install.packages("DMwR"); library(DMwR)
regnew <- SMOTE(CAN ~., data1, perc.over = 100, perc.under = 200)
table(regnew$CAN)
write.csv(regnew, "D:/SMOTEdata1.csv", row.names = F)
df1 <- read.csv("D:/SMOTEdata1.csv")

table(data2$CAN)
regnew <- SMOTE(CAN ~., data2, perc.over = 100, perc.under = 200)
table(regnew$CAN)
write.csv(regnew, "D:/SMOTEdata2.csv", row.names = F)
df2 <- read.csv("D:/SMOTEdata2.csv")


table(data3$CAN)
regnew <- SMOTE(CAN ~., data3, perc.over = 100, perc.under = 200)
table(regnew$CAN)
write.csv(regnew, "D:/SMOTEdata3.csv", row.names = F)
df3 <- read.csv("D:/SMOTEdata3.csv")
table(data4$CAN)
regnew <- SMOTE(CAN ~., data4, perc.over = 100, perc.under = 200)
table(regnew$CAN)
write.csv(regnew, "D:/SMOTEdata4.csv", row.names = F)
df4 <- read.csv("D:/SMOTEdata4.csv")

table(data5$CAN)
regnew <- SMOTE(CAN ~., data5, perc.over = 100, perc.under = 200)
table(regnew$CAN)
write.csv(regnew, "D:/SMOTEdata5.csv", row.names = F)
df5 <- read.csv("D:/SMOTEdata5.csv")

df<-df1 
#####################################################################

library(Hmisc)
matrix <- df %>%
  select(-Cancer)%>%
  as.matrix()

cor_matrix <- rcorr(matrix, type = "spearman")
cor_matrix$r
cor_matrix$P

rho <- cor_matrix$r
p_values <- cor_matrix$P

p_adjusted_fdr <- matrix(p.adjust(p_values, method = "fdr"), nrow = nrow(p_values))
rownames(p_adjusted_fdr) <- rownames(p_values)
colnames(p_adjusted_fdr) <- colnames(p_values)


signifiCancert_pairs <- which(p_adjusted_fdr < 0.05, arr.ind = TRUE)

for (i in 1:nrow(signifiCancert_pairs)) {
  row_idx <- signifiCancert_pairs[i, 1]
  col_idx <- signifiCancert_pairs[i, 2]
  var1 <- rownames(p_adjusted_fdr)[row_idx]
  var2 <- colnames(p_adjusted_fdr)[col_idx]
  
  if (row_idx < col_idx) {
    p_val <- p_adjusted_fdr[row_idx, col_idx]
    
    if (p_val < 0.001) {
      p_display <- "<0.001"
    } else {
      p_display <- format(round(p_val, 4), nsmall = 4)  
    }
    
    cat("显著相关对：", var1, "vs", var2, 
        "| Spearman ρ =", round(rho[row_idx, col_idx], 3), 
        "| 校正p值 =", p_display, "\n")
  }
}

library(corrplot)
dim(rho)
dim(p_adjusted_fdr)

corrplot(
  rho,
  p.mat = p_adjusted_fdr,
  type = "upper",          
  diag = FALSE,           
  insig = "blank",        
  sig.level = 0.05,
  method = "color",
  col = colorRampPalette(c("#1f98ff", "white", "#c01111"))(100),  
  tl.cex = 1,           
  tl.col = "black",       
  tl.srt = 45,            
  cl.cex = 1.5,           
  addCoef.col = "black",  
  number.cex = 0.8,       
  number.digits = 2,      
  mar = c(0.5, 0.5, 3.5, 0.5)    
)



library(car)
library(ggplot2)

model <- glm(
  Cancer ~ ., 
  data = df, 
  family = binomial(link = "logit")
)
vif_results <- vif(model)
print(vif_results)


vif_df <- as.data.frame(vif_results) %>% 
  mutate(Variable = rownames(.)) %>%  
  rename(GVIF = GVIF, Df = Df, Adjusted_VIF = "GVIF^(1/(2*Df))") %>% 
  dplyr::arrange(-Adjusted_VIF)  

print(vif_df)

ggplot(vif_df, aes(x = reorder(Variable, Adjusted_VIF), y = Adjusted_VIF, fill = Adjusted_VIF)) +
  geom_col() +
  geom_hline(yintercept = c(5, 10), linetype = "dashed", color = "pink") +
  scale_fill_gradient(low = "#1f98ff", high = "#d14") +
  labs(
    title = "Generalized Variance Inflation Factor (GVIF^{1/(2*Df)})",
    subtitle = "Diagnosis of Multicollinearity in overweight/obesity Prediction Model",
    x = "Predictor Variables",
    y = "Adjusted VIF Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(color = "#444444", size = 13, hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"), 
    panel.grid.major.y = element_blank(),
    
    axis.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, face = "bold", color = "black"),
    
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

#########################################################
library(mlr3)
library(dplyr)
library(precrec)
library(praznik)
library(mlr3verse)
library(mlr3extralearners) # remotes::install_github("mlr-org/mlr3extralearners")
library(mlr3proba)
library(mlr3tuning)
library(mlr3filters)
library(shapviz)
library(caret)

train_indices <- createDataPartition(df$Cancer, p = 0.6, list = FALSE)
train_data <- df[train_indices, ]
temp_data <- df[-train_indices, ]

valid_indices <- createDataPartition(temp_data$Cancer, p = 0.5, list = FALSE)
validation_data <- temp_data[valid_indices, ]
test_data <- temp_data[-valid_indices, ]

train_task <- as_task_classif(train_data, target = "Cancer", positive = "1")
test_task <- as_task_classif(test_data, target = "Cancer", positive = "1")
validation_task <- as_task_classif(validation_data, target = "Cancer", positive = "1")

library(xgboost)
xgb_lrn <- lrn("classif.xgboost",
               predict_type = "prob",
               nrounds = 500,  
               colsample_bytree = to_tune(lower = 0.1, upper = 1),  
               subsample =to_tune(lower = 0.1, upper = 1),  
               eta=to_tune(lower = 0.1, upper = 1),   
               max_depth = to_tune(lower = 1, upper = 10))

instance <- mlr3tuning::tune(tuner = tnr("grid_search"),
                             task = train_task,          
                             learner = xgb_lrn,
                             resampling = rsmp("cv", folds = 5), 
                             terminator = trm("evals", n_evals = 30, k = 0),
                             measures = msr("classif.auc"))

instance$result_learner_param_vals 
xgb_lrn$param_set$values <- instance$result_learner_param_vals
model_xgboost<-xgb_lrn$train(train_task) 


lgb_lrn <- lrn("classif.lightgbm",
               predict_type = "prob",
               learning_rate = to_tune(lower = 0.01, upper = 1), 
               bagging_fraction = to_tune(lower = 0.1, upper = 1),
               min_data_in_leaf=to_tune(lower = 10, upper = 50),
               max_depth = to_tune(lower = 1, upper = 10))
instance <- tune(tuner = tnr("grid_search"),
                 task = train_task,      
                 learner = lgb_lrn,
                 resampling = rsmp("cv", folds = 5), 
                 terminator = trm("evals", n_evals = 30, k = 0),
                 measures = msr("classif.auc"))

instance$result_learner_param_vals
lgb_lrn$param_set$values <- instance$result_learner_param_vals 
model_lightgbm<-lgb_lrn$train(train_task)

rf_lrn <- lrn("classif.randomForest",
              predict_type = "prob", 
              ntree = 500,           
              nodesize =to_tune(lower = 1, upper = 10), 
              importance="accuracy",
              mtry =to_tune(lower = 1, upper = 10))   


instance <- tune(tuner = tnr("grid_search"),
                 task = train_task,          
                 learner = rf_lrn,
                 resampling = rsmp("cv", folds = 5), 
                 terminator = trm("evals", n_evals = 30, k = 0),
                 measures = msr("classif.auc"))   

instance$result_learner_param_vals 

rf_lrn$param_set$values <- instance$result_learner_param_vals
model_rf<-rf_lrn$train(train_task)

knn_lrn <- lrn("classif.kknn",
               predict_type = "prob",            
               k = to_tune(lower = 1, upper = 30),     
               distance = to_tune(lower = 1, upper = 5),  
               kernel = to_tune(levels = c("rectangular", "gaussian"))) 

instance <- tune(
  tuner = tnr("grid_search"),
  task = train_task,
  learner = knn_lrn,  
  resampling = rsmp("cv", folds = 5),
  terminator = trm("evals", n_evals = 30),
  measures = msr("classif.auc")
)

instance$result_learner_param_vals 
knn_lrn$param_set$values <- instance$result_learner_param_val
model_knn<-knn_lrn$train(train_task) 


nnet_lrn <- lrn(
  "classif.nnet",
  predict_type = "prob",
  size = to_tune(lower = 3, upper = 10),
  decay = to_tune(lower = 0.0001, upper = 0.05)
  maxit = to_tune(lower = 50, upper = 200) 
)

instance <- tune(
  tuner = tnr("grid_search"), 
  task = train_task,
  learner = nnet_lrn,
  resampling = rsmp("cv", folds = 5), 
  terminator = trm("evals", n_evals = 30, k = 0), 
  measures = msr("classif.auc") 
)

instance$result_learner_param_vals
nnet_lrn$param_set$values <- instance$result_learner_param_vals 
model_nnet<-nnet_lrn$train(train_task)

##################################################
library(mlr3verse)
library(data.table)
library(pROC)
library(ggplot2)

model_colors <- c(
  "LightGBM" = "#1F77B4",
  "XGBoost" = "#FF7F0E", 
  "Random Forest" = "#2CA02C", 
  "KNN" = "#f20098",
  "NNET" ="yellow"
)

get_roc_data <- function(learner, task) {
  pred <- learner$predict(task)
  roc_obj <- roc(response = truth, predictor = probs, direction = "<", levels = c(0, 1))
  ci_sens <- ci.se(roc_obj, specificities = roc_obj$specificities)
  
  data.table(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    TPR_lower = ci_sens[, 1],
    TPR_upper = ci_sens[, 3],
    AUC = auc(roc_obj),
    # 修改后代码
    AUC_ci = sprintf("%.3f-%.3f", ci.auc(roc_obj)[1], ci.auc(roc_obj)[3]) # AUC的95% CI
  )
}


roc_dv <- rbindlist(list(
  get_roc_data(lgb_lrn, validation_task2)[, Model := "LightGBM"],
  get_roc_data(xgb_lrn, validation_task2)[, Model := "XGBoost"],
  get_roc_data(rf_lrn, validation_task2)[, Model := "Random Forest"],
  get_roc_data(knn_lrn, validation_task2)[, Model := "KNN"],
  get_roc_data(nnet_lrn, validation_task2)[, Model := "NNET"]
), use.names = TRUE)

ggplot(roc_dv, aes(x = FPR, y = TPR, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = TPR_lower, ymax = TPR_upper), alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 1.0, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(
    name = "Model (AUC [95% CI])",
    values = model_colors,
    # 修改后代码
    labels = ~ {
      auc_info <- unique(roc_dv[, .(Model, AUC, AUC_ci)])
      setkey(auc_info, Model)
      sprintf("%s: %.3f [%s]", auc_info$Model, auc_info$AUC, auc_info$AUC_ci)
    }
  ) +
  scale_fill_manual(values = model_colors, guide = "none") +  # 隐藏填充图例
  coord_equal(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "ROC Curve Comparison with 95% Confidence Intervals(Validation)",
    subtitle = "Ensemble Model Performance Benchmark",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = c(0.75, 0.2),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )

library(mlr3)
library(mlr3measures)

measures = msrs(c(
  "classif.auc",         
  "classif.acc",         
  "classif.precision",   
  "classif.recall",     
  "classif.specificity", 
  "classif.sensitivity",
  "classif.fpr",        
  "classif.fnr",         
  "classif.npv",         
  "classif.fbeta"       
))

metrics_lgb <- pred_lgb_v$score(measures)
metrics_xgb <- pred_xgb_v$score(measures)
metrics_rf <- pred_rf1_v$score(measures)

metrics_df <- rbind(
  data.frame(model = "LightGBM", as.list(metrics_lgb)),
  data.frame(model = "XGBoost", as.list(metrics_xgb)),
  data.frame(model = "Random Forest", as.list(metrics_rf))
)

colnames(metrics_df) <- c("model", "AUC", "Accuracy", "PPV", "Recall","Specificity","Sensitivity","FPR","FNR","NPV","F1-Score")
print(metrics_df)

library(ggplot2)
library(dplyr)
library(tidyr)
metric_order <- c("ROC AUC", "Accuracy", "PPV", "Recall","Specificity","FPR","FNR","NPV","F1-Score")

library(fmsb)
radar_df <- metrics_df[, -1]
colnames(radar_df) <- c("AUC", "Accuracy", "PPV", "Recall", "Specificity", "FPR", "FNR", "NPV", "F1")

radar_data <- rbind(
  rep(1, 9),   
  rep(0, 9), 
  radar_df
)
colors <- c("#1F77B4", "#FF7F0E", "#2CA02C")  

radarchart(
  radar_data,
  axistype = 1,            
  pcol = colors,           
  plwd = 2,               
  plty = 1,                
  cglcol = "grey",         
  cglty = 1,               
  cglwd = 0.8,             
  axislabcol = "grey",     
  vlcex = 0.8
)

legend(
  x = 1.2,          
  y = -0.3,         
  legend = metrics_df$model, 
  col = c("#1F77B4", "#FF7F0E", "#2CA02C"),
  lty = 1,
  lwd = 2,
  cex = 0.5,
  text.width = 0.6, 
  bty = "o",
  box.col = "gray",
  title = "Models",
  xpd = TRUE      
)
##########################################################

library(glmnet)
library(mice)
x <- train_data[, -which(names(train_data) == "Cancer")]
y <- train_data[, which(names(train_data) == "Cancer")]
fit <- glmnet(x, y, family = "binomial", alpha = 1)
print(fit) 
plot(fit, xvar = 'lambda', label = T)
x <- as.matrix(x)
cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = 'deviance')
print(cvfit)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
coef <- coef(cvfit, s = "lambda.1se")
coef <- coef@Dimnames[[1]][which(!coef == 0)]
coef[-1]
length(coef[-1])

library(Boruta)
BorutaOnIrisE <- Boruta(Cancer ~ ., data = train_data)
borutatable <- attStats(BorutaOnIrisE)
adjustdata <- function(data) {
  data <- cbind(rownames(data), data)
}
borutatable <- adjustdata(borutatable)
borutatable <- as.data.frame(borutatable)
names(borutatable) <- c("rowname", "meanimp", "medianimp", "minimp", "maximp", "normhits", "decision")
boruta_plot_df <- borutatable[order(borutatable$meanimp), ]
boruta_plot_df$rowname <- factor(boruta_plot_df$rowname, 
                                 levels = boruta_plot_df$rowname) 

decision_colors <- c("Confirmed" = "#f3993a",  
                     "Tentative" = "#108b96",  
                     "Rejected"  = "#E41A1C")  


boruta_plot_df$decision <- factor(boruta_plot_df$decision, 
                                  levels = c("Confirmed", "Tentative", "Rejected"))


confirmed_count <- sum(boruta_plot_df$decision == "Confirmed", na.rm = TRUE)
tentative_count <- sum(boruta_plot_df$decision == "Tentative", na.rm = TRUE)
rejected_count <- sum(boruta_plot_df$decision == "Rejected", na.rm = TRUE)


legend_labels <- c(
  paste0("Confirmed (", confirmed_count, ")"),
  paste0("Tentative (", tentative_count, ")"),
  paste0("Rejected (", rejected_count, ")")
)


ggplot(boruta_plot_df, aes(x = meanimp, y = rowname, fill = decision)) +
  geom_col() +
  scale_fill_manual(
    values = decision_colors,
    name = "Decision Status",
    breaks = c("Confirmed", "Tentative", "Rejected"),
    labels = legend_labels,  
    drop = FALSE,
    limits = c("Confirmed", "Tentative", "Rejected")
  ) +
  labs(
    title = "Boruta Feature Importance Analysis",
    subtitle = paste0("All ", confirmed_count, " variables are confirmed as important features"),
    x = "Mean Importance",
    y = "Variables"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    legend.position = "bottom"
  )

library(caret)
split <- createDataPartition(train_data$Cancer, p = 0.7, list = FALSE)
trainData <- train_data[split, ]
testData <- train_data[-split, ]
trainData$Cancer <- as.factor(trainData$Cancer)
control <- rfeControl(functions = rfFuncs, 
                      method = "cv",
                      number = 10)
sizes <- c(5,7,11,13,15,17,19,21) 
rfe_results <- rfe(x = trainData[, -ncol(trainData)], y = trainData$Cancer, sizes = sizes, rfeControl = control)
print(rfe_results)
plot(rfe_results, type = c("g", "o"))
selected_features <- predictors(rfe_results) 
print(selected_features)

results_df <- as.data.frame(rfe_results$results)
library(ggplot2)
library(ggthemes)
ggplot(results_df, aes(x = Variables, y = Accuracy, size = AccuracySD, color = Kappa)) +
  geom_point(alpha = 0.6) +
  scale_color_gradient(low = "#f3993a", high = "#108b96") + 
  scale_size(range = c(3, 15)) + 
  labs(title = "Recursive Feature Elimination Results", x = "Number of Features", y = "Accuracy", caption = "Bubble size indicates Accuracy SD; Color indicates Kappa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), plot.caption = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

library(dplyr)
lasso_vars  <- names(lasso_score)[lasso_score > 0]
boruta_vars <- borutatable$rowname[borutatable$decision == "Confirmed"]
rfe_vars    <- predictors(rfe_results)

all_vars <- colnames(train_data)[colnames(train_data) != "Cancer"]
total_n <- length(all_vars)

voting_df <- data.frame(Feature = all_vars) %>%
  mutate(
    Lasso_Hit  = as.numeric(Feature %in% lasso_vars),
    Boruta_Hit = as.numeric(Feature %in% boruta_vars),
    RFE_Hit    = as.numeric(Feature %in% rfe_vars),
    
    Rank_Lasso = ifelse(Feature %in% names(lasso_score), 
                        rank(-lasso_score, ties.method = "min")[Feature], 
                        total_n),
    
    Rank_Boruta = ifelse(Feature %in% names(boruta_score), 
                         rank(-boruta_score, ties.method = "min")[Feature], 
                         total_n),
    
    Rank_RFE = ifelse(Feature %in% names(rfe_score), 
                      rank(-rfe_score, ties.method = "min")[Feature], 
                      total_n)
  ) %>%
  mutate(
    Frequency = Lasso_Hit + Boruta_Hit + RFE_Hit,
    Rank_Sum  = Rank_Lasso + Rank_Boruta + Rank_RFE
  )

voting_df <- voting_df %>%
  arrange(desc(Frequency), Rank_Sum)

voting_df$Consensus_Rank <- 1:nrow(voting_df)

library(mlr3)
library(mlr3learners)

results_consensus <- data.frame()
for (i in 1:length(consensus_features)) {
  current_selection <- consensus_features[1:i]
  task_tmp <- train_task$clone()$select(current_selection)
  test_task_tmp <- test_task$clone()$select(current_selection)
  xgb_lrn$train(task_tmp)
  pred <- xgb_lrn$predict(test_task_tmp)
  scores <- pred$score(msrs(c("classif.acc", "classif.auc", "classif.fbeta")))

  results_consensus <- rbind(results_consensus, data.frame(
    num_features = i,
    accuracy = scores[["classif.acc"]],
    auroc = scores[["classif.auc"]],
    f1 = scores[["classif.fbeta"]]
  ))
}


results_long <- melt(results_consensus, id.vars = "num_features", 
                     variable.name = "metric", value.name = "score")
final_line_plot <- ggplot(results_long, aes(x = num_features, y = score, color = metric)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(aes(shape = metric), size = 3) +
  scale_color_jama(labels = c("Accuracy", "AUROC", "F1 Score")) +
  scale_shape_manual(values = c(17, 16, 15), labels = c("Accuracy", "AUROC", "F1 Score")) +
  scale_x_continuous(breaks = seq(1, max(results_long$num_features), by = 1)) +
  labs(
    title = "Model Performance via Incremental Feature Selection",
    subtitle = "Features ranked by consensus of LASSO, Boruta, and RFE",
    x = "Number of Top-ranked Features Included",
    y = "Performance Metric Value"
  ) 

###################################
final_features <- consensus_features[1:7]

train_task_final <- train_task$clone()$select(final_features)
test_task_final  <- test_task$clone()$select(final_features)
val_task_final   <- validation_task$clone()$select(final_features)

model_xgboost_final <- xgb_lrn$train(train_task_final)

library(shapviz)
library(ggplot2)
X_train_mat <- as.matrix(train_data[, final_features])
native_xgb_final <- model_xgboost_final$model
shap_train <- shapviz(native_xgb_final, X_pred = X_train_mat)
sv_importance(shap_train, kind = "beeswarm", max_display = 10) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 14, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
sv_importance(shap_train, kind = "bar", max_display = 10) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold")
  ) 


library(lime)
library(MASS)
library(caret)
library(e1071)
library(lime)
library(dplyr)
lime_train_data <- as.data.frame(train_data[, final_features])
lime_test_data  <- as.data.frame(test_data[, final_features])
native_xgb <- model_xgboost_final$model

predict_model.xgb.Booster <- function(x, newdata, type, ...) {
  m <- xgb.DMatrix(as.matrix(newdata))
  res <- predict(x, m)
  return(data.frame(Cancer = res, Healthy = 1 - res))
}

model_type.xgb.Booster <- function(x, ...) {
  return("classification")
}
explainer <- lime(
  x = lime_train_data, 
  model = native_xgb,
  bin_continuous = TRUE,   
  n_bins = 4
)

explanation <- explain(
  x = lime_test_data[1:6, ], 
  explainer = explainer, 
  n_labels = 1, 
  n_features = 7
)
plot_features(explanation) +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    title = element_text(size = 14, face = "bold")
  ) +
  labs(title = "LIME Individual Feature Importance",
       subtitle = "Local explanations for the first 3 test samples")

