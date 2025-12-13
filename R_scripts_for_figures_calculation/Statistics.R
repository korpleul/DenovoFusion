# In this file, we apply such statistics methods on check results from 4 tools:
# HMFtools, FACTERA, Genefuse, DenovoFusion.
# normal evaluation: sensitivity, specificify, accuracy F1
# paired comparison: McNemar test, Cohan Kappa test with Benjamini–Hochberg false discovery rate (FDR) correction
# all groups in logistic regression model



library(tidyr)
library(dplyr)
library(clipr)
library(broom)
library(ggplot2)
library(dplyr)
library(psych)

# Reading Data
df <- read.csv("fusion_results.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert to long format (each line = one sample and one tool)
df_long <- pivot_longer(df, 
                        cols = -Samples,
                        names_to = "Tool", 
                        values_to = "Category")

head(df_long)

# Defining statistical functions
compute_metrics <- function(pred_vector) {
  TP <- sum(pred_vector == "TP")
  TN <- sum(pred_vector == "TN")
  FP <- sum(pred_vector == "FP")
  FN <- sum(pred_vector == "FN")
  
  Sensitivity <- TP / (TP + FN)
  Specificity <- TN / (TN + FP)
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  Precision <- TP / (TP + FP)
  F1 <- 2 * (Precision * Sensitivity) / (Precision + Sensitivity)
  
  return(data.frame(TP, TN, FP, FN, Sensitivity, Specificity, Accuracy, F1))
}
# Calculate for each tool
tools <- colnames(df)[-1]  # Exclude the Sample column
metrics_list <- lapply(tools, function(tool) compute_metrics(df[[tool]]))
names(metrics_list) <- tools

# Merge into one table
metrics_df <- do.call(rbind, metrics_list)
metrics_df <- cbind(Tool = rownames(metrics_df), metrics_df)
rownames(metrics_df) <- NULL

# View Results
print(metrics_df)
write_clip(metrics_df)




# 1111
# McNemar’s test, p = 4.19e-11
# Cohen's Kappa between HMFtools and FACTERA was 0.161 (p = 0.428e-09), indicating low agreement
# Extract the results of both tools
tool_A <- df_long %>% filter(Tool == "HMFtools") %>% arrange(Samples)
tool_B <- df_long %>% filter(Tool == "FACTERA") %>% arrange(Samples)


# McNemar test
# Defining criteria for correct predictions (TP and TN)
correct_A <- tool_A$Category %in% c("TP", "TN")
correct_B <- tool_B$Category %in% c("TP", "TN")

# Constructing a 2x2 cross table
confusion_table <- table(correct_A, correct_B)
print(confusion_table)

# Use McNemar's test
mcnemar.test(confusion_table)



# Cohen's Kappa test
# Make sure the number of rows is the same
stopifnot(nrow(tool_A) == nrow(tool_B))

# Construct a data frame containing the classification results of the two tools
kappa_input <- data.frame(tool_A = tool_A$Category, tool_B = tool_B$Category)

# Cohen's Kappa
kappa_result <- cohen.kappa(kappa_input)

print(kappa_result)

kappa_value <- kappa_result$kappa
lower_ci <- kappa_result$confid[1,1]
upper_ci <- kappa_result$confid[1,3]
se <- sqrt(kappa_result$var.kappa)
z <- kappa_value / se
p_value <- 2 * (1 - pnorm(abs(z)))


kappa_df <- data.frame(
  ToolPair = "Tool A vs Tool B",
  Kappa = round(kappa_value, 3),
  LowerCI = round(lower_ci, 3),
  UpperCI = round(upper_ci, 3),
  P_value = signif(p_value, 3),
  Signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
)

print(kappa_df)





# 2222
# McNemar’s test, p = 0.1547
# Cohen's Kappa between HMFtools and Genefuse was 0.340 (p = 4.528e-07), indicating low agreement.
# Extract the results of both tools
tool_A <- df_long %>% filter(Tool == "HMFtools") %>% arrange(Samples)
tool_B <- df_long %>% filter(Tool == "Genefuse") %>% arrange(Samples)

# Defining criteria for correct predictions (TP and TN)
correct_A <- tool_A$Category %in% c("TP", "TN")
correct_B <- tool_B$Category %in% c("TP", "TN")

# Constructing a 2x2 cross table
confusion_table <- table(correct_A, correct_B)
print(confusion_table)

# Use McNemar's test
mcnemar.test(confusion_table)

# Cohen's Kappa test
# Make sure the number of rows is the same
stopifnot(nrow(tool_A) == nrow(tool_B))

# Construct a data frame containing the classification results of the two tools
kappa_input <- data.frame(tool_A = tool_A$Category, tool_B = tool_B$Category)

# Cohen's Kappa
kappa_result <- cohen.kappa(kappa_input)


print(kappa_result)


kappa_value <- kappa_result$kappa
lower_ci <- kappa_result$confid[1,1]
upper_ci <- kappa_result$confid[1,3]
se <- sqrt(kappa_result$var.kappa)
z <- kappa_value / se
p_value <- 2 * (1 - pnorm(abs(z)))


kappa_df <- data.frame(
  ToolPair = "Tool A vs Tool B",
  Kappa = round(kappa_value, 3),
  LowerCI = round(lower_ci, 3),
  UpperCI = round(upper_ci, 3),
  P_value = signif(p_value, 3),
  Signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
)

print(kappa_df)


# 3333
# McNemar's chi-squared = 15.848, df = 1, p-value = 6.865e-05
# Cohen's Kappa between HMFtools and DenovoFusion was 0.337 (p = 2.831e-09), indicating low agreement.
# Extract the results of both tools
tool_A <- df_long %>% filter(Tool == "HMFtools") %>% arrange(Samples)
tool_B <- df_long %>% filter(Tool == "DenovoFusion") %>% arrange(Samples)

# Defining criteria for correct predictions (TP and TN)
correct_A <- tool_A$Category %in% c("TP", "TN")
correct_B <- tool_B$Category %in% c("TP", "TN")

# Constructing a 2x2 cross table
confusion_table <- table(correct_A, correct_B)
print(confusion_table)

# Use McNemar's test
mcnemar.test(confusion_table)

# Cohen's Kappa test
# Make sure the number of rows is the same
stopifnot(nrow(tool_A) == nrow(tool_B))

# Construct a data frame containing the classification results of the two tools
kappa_input <- data.frame(tool_A = tool_A$Category, tool_B = tool_B$Category)

# Cohen's Kappa
kappa_result <- cohen.kappa(kappa_input)


print(kappa_result)

kappa_value <- kappa_result$kappa
lower_ci <- kappa_result$confid[1,1]
upper_ci <- kappa_result$confid[1,3]
se <- sqrt(kappa_result$var.kappa)
z <- kappa_value / se
p_value <- 2 * (1 - pnorm(abs(z)))


kappa_df <- data.frame(
  ToolPair = "Tool A vs Tool B",
  Kappa = round(kappa_value, 3),
  LowerCI = round(lower_ci, 3),
  UpperCI = round(upper_ci, 3),
  P_value = signif(p_value, 3),
  Signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
)

print(kappa_df)

# 4444
# McNemar's chi-squared = 35.2, df = 1, p-value = 2.975e-09
# Cohen's Kappa between FACTERA and Genefuse was 0.221 (p = 7.717e-11), indicating low agreement.
# Extract the results of both tools
tool_A <- df_long %>% filter(Tool == "FACTERA") %>% arrange(Samples)
tool_B <- df_long %>% filter(Tool == "Genefuse") %>% arrange(Samples)

# Defining criteria for correct predictions (TP and TN)
correct_A <- tool_A$Category %in% c("TP", "TN")
correct_B <- tool_B$Category %in% c("TP", "TN")

# Constructing a 2x2 cross table
confusion_table <- table(correct_A, correct_B)
print(confusion_table)

# Use McNemar's test
mcnemar.test(confusion_table)

# Cohen's Kappa test
# Make sure the number of rows is the same
stopifnot(nrow(tool_A) == nrow(tool_B))

# Construct a data frame containing the classification results of the two tools
kappa_input <- data.frame(tool_A = tool_A$Category, tool_B = tool_B$Category)

# Cohen's Kappa
kappa_result <- cohen.kappa(kappa_input)


print(kappa_result)


kappa_value <- kappa_result$kappa
lower_ci <- kappa_result$confid[1,1]
upper_ci <- kappa_result$confid[1,3]
se <- sqrt(kappa_result$var.kappa)
z <- kappa_value / se
p_value <- 2 * (1 - pnorm(abs(z)))


kappa_df <- data.frame(
  ToolPair = "Tool A vs Tool B",
  Kappa = round(kappa_value, 3),
  LowerCI = round(lower_ci, 3),
  UpperCI = round(upper_ci, 3),
  P_value = signif(p_value, 3),
  Signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
)

print(kappa_df)

# 5555
# McNemar's chi-squared = 25.037, df = 1, p-value = 5.624e-07
# Cohen's Kappa between FACTERA and DenovoFusion was 0.467 (p = 2.22e-16), indicating moderate agreement.
# Extract the results of both tools
tool_A <- df_long %>% filter(Tool == "FACTERA") %>% arrange(Samples)
tool_B <- df_long %>% filter(Tool == "DenovoFusion") %>% arrange(Samples)

# Defining criteria for correct predictions (TP and TN)
correct_A <- tool_A$Category %in% c("TP", "TN")
correct_B <- tool_B$Category %in% c("TP", "TN")

# Constructing a 2x2 cross table
confusion_table <- table(correct_A, correct_B)
print(confusion_table)

# Use McNemar's test
mcnemar.test(confusion_table)

# Cohen's Kappa test
# Make sure the number of rows is the same
stopifnot(nrow(tool_A) == nrow(tool_B))

# Construct a data frame containing the classification results of the two tools
kappa_input <- data.frame(tool_A = tool_A$Category, tool_B = tool_B$Category)

# Cohen's Kappa
kappa_result <- cohen.kappa(kappa_input)


print(kappa_result)


kappa_value <- kappa_result$kappa
lower_ci <- kappa_result$confid[1,1]
upper_ci <- kappa_result$confid[1,3]
se <- sqrt(kappa_result$var.kappa)
z <- kappa_value / se
p_value <- 2 * (1 - pnorm(abs(z)))


kappa_df <- data.frame(
  ToolPair = "Tool A vs Tool B",
  Kappa = round(kappa_value, 3),
  LowerCI = round(lower_ci, 3),
  UpperCI = round(upper_ci, 3),
  P_value = signif(p_value, 3),
  Signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
)

print(kappa_df)

# 6666
# McNemar's chi-squared = 7.225, df = 1, p-value = 0.00719
# Cohen's Kappa between Genefuse and DenovoFusion was 0.394 (p = 1.084e-09), indicating low agreement.
# Extract the results of both tools
tool_A <- df_long %>% filter(Tool == "Genefuse") %>% arrange(Samples)
tool_B <- df_long %>% filter(Tool == "DenovoFusion") %>% arrange(Samples)

# Defining criteria for correct predictions (TP and TN)
correct_A <- tool_A$Category %in% c("TP", "TN")
correct_B <- tool_B$Category %in% c("TP", "TN")

# Constructing a 2x2 cross table
confusion_table <- table(correct_A, correct_B)
print(confusion_table)

# Use McNemar's test
mcnemar.test(confusion_table)

# Cohen's Kappa test
# Make sure the number of rows is the same
stopifnot(nrow(tool_A) == nrow(tool_B))

# Construct a data frame containing the classification results of the two tools
kappa_input <- data.frame(tool_A = tool_A$Category, tool_B = tool_B$Category)

# Cohen's Kappa
kappa_result <- cohen.kappa(kappa_input)


print(kappa_result)


kappa_value <- kappa_result$kappa
lower_ci <- kappa_result$confid[1,1]
upper_ci <- kappa_result$confid[1,3]
se <- sqrt(kappa_result$var.kappa)
z <- kappa_value / se
p_value <- 2 * (1 - pnorm(abs(z)))


kappa_df <- data.frame(
  ToolPair = "Tool A vs Tool B",
  Kappa = round(kappa_value, 3),
  LowerCI = round(lower_ci, 3),
  UpperCI = round(upper_ci, 3),
  P_value = signif(p_value, 3),
  Signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
)

print(kappa_df)





### Figures
base_theme <- function(base_size = 16) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
   
      text = element_text(face = "bold", color = "black"),
      
      
      axis.text = element_text(size = rel(1.2)),  
      axis.title = element_text(size = rel(1.2)),
      axis.text.x = element_text(margin = margin(t = 1)),  
      axis.text.y = element_text(margin = margin(r = 1)),
      
      legend.text = element_text(size = rel(1.0)),
      legend.title = element_text(size = rel(1)),
      
      
      plot.title   = element_text(size = rel(1.3), hjust = 0.5),
      plot.caption = element_text(size = rel(1)),
      
      strip.text = element_text(size = rel(1.2))
    )
}


# logistic regression model and ouTPut figure
df_long <- df_long %>%
  mutate(Correct = ifelse(Category %in% c("TP", "TN"), 1, 0))


# First set the factor order and specify the reference group
tools_order <- c("DenovoFusion", "HMFtools", "FACTERA", "Genefuse")
df_long$Tool <- factor(df_long$Tool, levels = tools_order)

df_long$Tool <- relevel(df_long$Tool, ref = "DenovoFusion")  

# Refit the model
model <- glm(Correct ~ Tool, data = df_long, family = "binomial")

summary(model)

# Calculate OR and confidence interval
exp(cbind(OR = coef(model), confint(model)))


# Calculate OR and confidence interval
model_df <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
model_df <- model_df %>% filter(term != "(Intercept)")

model_df$term <- recode(model_df$term,
                        "ToolFACTERA" = "FACTERA",
                        "ToolGenefuse" = "Genefuse",
                        "ToolHMFtools" = "HMFtools")

# Calculate the prediction accuracy of each tool
accuracy_df <- df_long %>%
  group_by(Tool) %>%
  summarise(Accuracy = mean(Correct) * 100) %>%
  arrange(desc(Accuracy))

# Merge accuracy into model_df
model_df <- left_join(model_df, accuracy_df, by = c("term" = "Tool"))


model_df <- model_df %>%
  mutate(signif = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE           ~ " "
  ))

# Figure
ggplot(model_df, aes(x = reorder(term, estimate), y = estimate)) +
  geom_bar(stat = "identity", fill = "#666666", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#7f7f7f") +
  geom_text(aes(label = sprintf("%.1f%%", Accuracy), y = conf.high + 0.15), color = "black", size = 6, family = "Arial", fontface = "bold") +  
  geom_text(aes(label = signif, y = conf.high + 0.4), color = "#d62728", size = 6, family = "Arial", fontface = "bold", size = 11/.pt) +                         
  labs(title = "Logistic Regression: Tools comparison)",
       x = NULL, y = "Odds Ratio (vs. DenovoFusion)") +
  base_theme(base_size = 18) +
  ylim(0, max(model_df$conf.high)*1.3)  


# McNemar bar chart
mcnemar_df <- data.frame(
  ToolPair = c("H/F", "H/G", "H/D", 
               "F/G", "F/D", "G/D"),
  DisagreeA = c(61, 25, 36, 5, 0, 29),   
  DisagreeB = c(6, 15, 9, 50, 28, 12),   
  Pvalue = c(0.0001, 0.11, 0.0001, 0.0001, 0.0001, 0.013),
  Signif = c("***", " ", "***", "***", "***", "*")
)

# FDR (Benjamini–Hochberg) correction
mcnemar_df$Pvalue_FDR <- p.adjust(mcnemar_df$Pvalue, method = "BH")

# Regenerate significance markers using FDR-corrected p-values
mcnemar_df$Signif <- case_when(
  mcnemar_df$Pvalue_FDR < 0.001 ~ "***",
  mcnemar_df$Pvalue_FDR < 0.01  ~ "**",
  mcnemar_df$Pvalue_FDR < 0.05  ~ "*",
  TRUE ~ " "
)

library(tidyr)
mcnemar_long <- pivot_longer(mcnemar_df, cols = c("DisagreeA", "DisagreeB"),
                             names_to = "Direction", values_to = "Count")
mcnemar_long$ToolPair <- factor(mcnemar_long$ToolPair, levels = mcnemar_df$ToolPair)


ggplot(mcnemar_long, aes(x = ToolPair, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"),
                    labels = c("A correct / B wrong", "A wrong / B correct")) +
  geom_text(data = mcnemar_df,
            aes(x = ToolPair, y = pmax(DisagreeA, DisagreeB) + 2,
                label = Signif),
            color = "#d62728", size = 6,
            size = 14/.pt,  
            family = "Arial",  
            fontface = "bold",  
            inherit.aes = FALSE) +
  labs(title = "McNemar's Test Disagreement Counts",
       x = NULL, y = "Count", fill = "Disagreement Type") +
  base_theme(base_size = 18) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
  legend.position = c(0.78, 0.92))  #  Adjust legend position





# Cohan kappa test dot plot + error bars

kappa_df <- data.frame(
  ToolPair = c("H/F", "H/G", "H/D", 
               "F/G", "F/D", "G/D"),
  Kappa = c(0.172, 0.356, 0.36, 0.221, 0.456, 0.377),
  LowerCI = c(0.0978, 0.2040, 0.234, 0.129, 0.316, 0.233),
  UpperCI = c(0.247, 0.508, 0.489, 0.312, 0.596, 0.521),
  Signif = c("***", "***", "***", "***", "***", "***")
)
kappa_df$ToolPair <- factor(kappa_df$ToolPair, levels = kappa_df$ToolPair)

ggplot(kappa_df, aes(x = ToolPair, y = Kappa)) +
  geom_point(size = 4, color = "#666666") +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.15, color = "#333333") +
  geom_text(aes(label = Signif, y = UpperCI + 0.02), color = "#d62728",  
            size = 14/.pt,         
            family = "Arial",      
            fontface = "bold") +
  coord_cartesian(ylim = c(min(kappa_df$LowerCI) - 0.05, 1)) +
  labs(title = "Cohen's Kappa Test of Tool Pairs", y = "Kappa Value", x = NULL) +
  base_theme(base_size = 18) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )


# bar chart for running time and memory usage


perf_df <- data.frame(
  Tool = c("HMFtools", "FACTERA", "Genefuse", "DenovoFusion"),
  Time_sec = c(5.5, 25.4, 1.1, 16.1),       
  MaxMemory_GB = c(32.2, 16.5, 12, 75)               
)


perf_df$Tool <- factor(perf_df$Tool, levels = perf_df$Tool)

nature_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

p1 <- ggplot(perf_df, aes(x = Tool, y = Time_sec, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = nature_colors) +
  scale_x_discrete(expand = c(0.05, 0)) + 
  labs(title = "Runtime of Fusion Detection Tools", y = "Time (hours)", x = NULL) +
  base_theme(base_size = 18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

p1

p2 <- ggplot(perf_df, aes(x = Tool, y = MaxMemory_GB, fill = Tool)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = nature_colors) +
  scale_x_discrete(expand = c(0.05, 0)) + 
  labs(title = "Peak Memory Usage", y = "Max Memory (GB)", x = NULL) +
  base_theme(base_size = 18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )


p2



# speedup and efficiency figure
# input the data
df <- data.frame(
  Tool = rep(c("HMFtools", "Genefuse", "FACETERA", "DenovoFusion"), each = 5),
  Threads = factor(rep(c(1, 2, 4, 8, 16), times = 4), levels = c(1, 2, 4, 8, 16)),
  Efficiency = c(
    1.00, 0.6, 0.52, 0.34, 0.24,   # HMFtools
    1.00, 0.71, 0.42, 0.22, 0.12,   # Genefuse
    1.00, 0.88, 0.65, 0.55, 0.42,   # FACETERA
    1.00, 0.85, 0.70, 0.55, 0.45    # DenovoFusion
  )
)


p3 <- ggplot(df, aes(x = Threads, y = Efficiency, color = Tool, group = Tool)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_y_reverse(limits = c(1, 0.1), breaks = seq(0, 1, 0.1)) +
  scale_color_manual(
    values = nature_colors,
    breaks = rev(levels(factor(df$Tool)))  # This reverses the legend order
  ) +
  labs(title = "Parallel Efficiency vs Number of Threads",
       x = "Number of Threads",
       y = "Efficiency (reversed)") +
  base_theme(base_size = 18) +
  theme(legend.position = c(0.05, 0.95),  # top-left inside plot (x, y as fraction)
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16))

p3
