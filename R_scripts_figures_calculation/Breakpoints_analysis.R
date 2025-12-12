# This script created the figures which showed the distribution of breakpoints in our samples
# See Figure5A-C

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(readr)
library(openxlsx)
library(ggridges)


# set a size for all figures
base_theme <- function(base_size = 16) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      
      text = element_text(face = "bold", color = "black"),
      
      
      axis.text = element_text(size = rel(1)),  
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

# input the data from different sheets
breakpoints_list <- list(
  "EWSR1-FLI1" = read_excel("breakpoints.xlsx", sheet = "EWSR1-FLI1"),
  "EWSR1-ERG"  = read_excel("breakpoints.xlsx", sheet = "EWSR1-ERG"),
  "EWSR1-ETV1" = read_excel("breakpoints.xlsx", sheet = "EWSR1-ETV1")
)

# Define functions to calculate statistics and normality tests
analyze_breakpoints <- function(df, fusion_name) {
  result <- data.frame(
    Fusion = character(),
    Position = character(),
    Mean = numeric(),
    SD = numeric(),
    Min = numeric(),
    Max = numeric(),
    Shapiro_W_pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (nrow(df) == 0) {
    warning(paste("Warning: Data frame for", fusion_name, "is empty. Skipping."))
    return(result)
  }
  
  for (pos_col in c("pos1", "pos2")) {
    vals <- df[[pos_col]]
    if (length(vals) < 3) {
      warning(paste("Warning: Not enough data points for", fusion_name, pos_col, "to perform Shapiro-Wilk test."))
      next
    }
    
    shapiro_res <- tryCatch(
      shapiro.test(vals),
      error = function(e) {
        warning(paste("Shapiro-Wilk test error for", fusion_name, pos_col, ":", e$message))
        return(list(p.value = NA))
      }
    )
    
    result <- rbind(result, data.frame(
      Fusion = fusion_name,
      Position = pos_col,
      Mean = mean(vals),
      SD = sd(vals),
      Min = min(vals),
      Max = max(vals),
      Shapiro_W_pvalue = shapiro_res$p.value,
      stringsAsFactors = FALSE
    ))
  }
  return(result)
}


# Calculate statistics for each fusion gene
all_results <- lapply(names(breakpoints_list), function(fn) {
  analyze_breakpoints(breakpoints_list[[fn]], fn)
}) %>% bind_rows()

print(all_results)

# Excel
write.xlsx(all_results, "breakpoints_statistics.xlsx")

# consistency statistics table
# Function to compute consistency statistics per fusion type
compute_consistency_summary <- function(df, fusion_name, max_diff = 50) {
  df_grouped <- df %>%
    group_by(Samples, Fusion_name) %>%
    summarise(
      pos1_range = max(pos1) - min(pos1),
      pos2_range = max(pos2) - min(pos2),
      n_methods = n(),
      .groups = "drop"
    ) %>%
    mutate(
      pos1_consistent = pos1_range <= max_diff,
      pos2_consistent = pos2_range <= max_diff,
      both_consistent = pos1_consistent & pos2_consistent
    )
  
  total_samples <- n_distinct(df$Samples)
  eligible_samples <- sum(df_grouped$n_methods > 1)
  
  stats <- df_grouped %>% filter(n_methods > 1)
  
  summary_row <- tibble(
    Fusion = fusion_name,
    Total_samples = total_samples,
    Eligible_samples = eligible_samples,
    Pos1_consistent_samples = sum(stats$pos1_consistent),
    Pos1_consistency_rate = round(mean(stats$pos1_consistent) * 100, 2),
    Pos2_consistent_samples = sum(stats$pos2_consistent),
    Pos2_consistency_rate = round(mean(stats$pos2_consistent) * 100, 2),
    Both_consistent_samples = sum(stats$both_consistent),
    Both_consistency_rate = round(mean(stats$both_consistent) * 100, 2)
  )
  
  return(summary_row)
}

# Apply to each fusion sheet
consistency_summary_table <- bind_rows(
  compute_consistency_summary(breakpoints_list$`EWSR1-FLI1`, "EWSR1-FLI1"),
  compute_consistency_summary(breakpoints_list$`EWSR1-ERG`, "EWSR1-ERG"),
  compute_consistency_summary(breakpoints_list$`EWSR1-ETV1`, "EWSR1-ETV1")
)

# Write output
write_csv(consistency_summary_table, "fusion_consistency_summary.csv")



# densityplot for breakpoints deviation
# Read EWSR1-FLI1 data
df <- read_excel("breakpoints.xlsx", sheet = "EWSR1-FLI1")

# Make sure the chr column is of type character (prefix "chr" to match the standard chromosome format)
df <- df %>%
  mutate(
    chr1 = paste0("chr", chr1),
    chr2 = paste0("chr", chr2)
  )

df$Methods <- factor(df$Methods, levels = c("HMFtools", "Genefuse", "DenovoFusion"))

# density（5'pos1）
ggplot(df, aes(x = pos1, fill = Methods)) +
  geom_histogram(binwidth = 100, color = "white", alpha = 0.8, position = "identity") +
  facet_wrap(~ Methods, ncol = 1, scales = "free_y") +
  labs(title = "5′ Breakpoint Histogram (pos1)", x = "Genomic Position (chr22)", y = "Count") +
  base_theme(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

# density （3'pos2）
ggplot(df, aes(x = pos2, fill = Methods)) +
  geom_histogram(binwidth = 600, color = "white", alpha = 0.8, position = "identity") +
  facet_wrap(~ Methods, ncol = 1, scales = "free_y") +
  labs(title = "3′ Breakpoint Histogram (pos2)", x = "Genomic Position (chr11)", y = "Count") +
  base_theme(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position = "none")


# Figures
methods_order <- unique(df$Methods)
df$Methods <- factor(df$Methods, levels = methods_order)

p <- ggplot(df, aes(x = pos1, y = pos2)) +
  geom_point(aes(color = Methods, alpha = compo1), size = 4) +
  coord_fixed(ratio = 0.2) +
  scale_x_continuous(
    name = "5′ Breakpoint (EWSR1 - chr22)",
    limits = c(29681000, 29690000),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    name = "3′ Breakpoint (FLI1 - chr11)",
    limits = c(128640000, 128686000),
    expand = expansion(mult = 0.02)
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_alpha_manual(values = c("exon" = 1, "intron" = 1), guide = "none") +
  labs(title = "Breakpoint Coordinate Plot for EWSR1-FLI1") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


ggMarginal(p, type = "density", fill = "grey", alpha = 0.5)


