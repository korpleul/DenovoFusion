# This file created the heatmap of the detected ETS-family related fusion genes from different methods
# see Figure3A

library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)

# load the data
df <- read_excel("fusion_results.xlsx", sheet = "Sheet3")

# convert to long dataframe
df_long <- pivot_longer(df, 
                        cols = -Samples,
                        names_to = "Tool", 
                        values_to = "Category")

# Four categories: EF(EWSR1-FLI1), EE(EWSR1-ERG), EV(EWSR1-ETV1), /(NO SUCH FUSION)
# Others are temporarily classified as "Other" for easy expansion
df_long <- df_long %>%
  mutate(Status = case_when(
    Category == "EWSR1-FLI1" ~ "EWSR1-FLI1",
    Category == "EWSR1-ERG" ~ "EWSR1-ERG",
    Category == "EWSR1-ETV1" ~ "EWSR1-ETV1",
    Category == "/" ~ "Not detected",
    TRUE ~ "Other"
  ))

# Set the factor order of Status
df_long <- df_long %>%
  mutate(Status = factor(Status, levels = c("EWSR1-FLI1", "EWSR1-ERG", "EWSR1-ETV1", "Other", "Not detected")))


# 4 catgories
category_colors <- c(
  "EWSR1-FLI1" = "#1f77b4",   
  "EWSR1-ERG" = "#ff7f0e",   
  "EWSR1-ETV1" = "#2ca02c",   
  "Not detected" = "#7f7f7f",   
  "Other" = "#F5F5F5" 
)


# To draw finer bars, convert samples and tools into factors and sort them
df_long$Samples <- factor(df_long$Samples, levels = unique(df_long$Samples))
df_long$Tool <- factor(df_long$Tool, levels = rev(unique(df_long$Tool)))  # Reverse the Y axis order for a more beautiful look

# Drawing
p <- ggplot(df_long, aes(x = Samples, y = Tool, fill = Status)) +
  geom_tile(color = "gray90", linewidth = 0.15) +
  scale_fill_manual(values = category_colors) +
  labs(title = "TP/TN Highlighted Fusion Detection (Comparison of Fusion Detection Tools)", 
       x = "", 
       y = " ") +
  theme_minimal(base_size = 24) +
  theme(
    text = element_text(face = "bold", color = "#000000", family = "Arial"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14,  margin = margin(t = 1)),
    axis.text.y = element_text(size = 24, hjust = 0),
    legend.text = element_text(size = 16),  
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 2.5)

# Display graphics
print(p)

ggsave("3A.svg", 
       plot = p,               
       width = 21,             
       height = 4,            
       dpi = 300,               )    
dev.off()


library(ggplot2)
library(tidyr)
library(dplyr)

# Reading Data
df <- read.csv("fusion_results.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert to long format
df_long <- pivot_longer(df, 
                        cols = -Samples,
                        names_to = "Tool", 
                        values_to = "Category")

# Four categories: TP, TN, FP, FN
# Your data corresponds to these four types of strings, make sure they match exactly
# Others are temporarily classified as "Other" for easy expansion
df_long <- df_long %>%
  mutate(Status = case_when(
    Category == "TP" ~ "TP",
    Category == "TN" ~ "TN",
    Category == "FP" ~ "FP",
    Category == "FN" ~ "FN",
    TRUE ~ "Other"
  ))

# Four categories of color matching
category_colors <- c(
  "TP" = "#08306B",  
  "TN" = "#DEEBF7",   
  "FP" = "#99000D",   
  "FN" = "#FCBBA1",   
  "Other" = "#F5F5F5" 
)


# To draw finer bars, convert samples and tools into factors and sort them
df_long$Samples <- factor(df_long$Samples, levels = unique(df_long$Samples))
df_long$Tool <- factor(df_long$Tool, levels = rev(unique(df_long$Tool)))  # Reverse the Y axis order for a more beautiful look

# Drawing
p <- ggplot(df_long, aes(x = Samples, y = Tool, fill = Status)) +
  geom_tile(color = "gray90", linewidth = 0.15) +
  scale_fill_manual(values = category_colors) +
  labs(title = "", x = "TP/TN Highlighted Fusion Detection(Comparison of Fusion Detection Tools) ", y = " ") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "left",
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +
  coord_fixed(ratio = 2.5)  

# Display graphics
print(p)

ggsave("Fusion_tools_heatmap.svg", 
       plot = p,               
       width = 16,              
       height = 4,            
       dpi = 96,               
       bg = "transparent")     
dev.off()
