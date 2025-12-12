# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite) # Essential for saving SVG files

# ==============================================================================
# 1. Prepare Data
# ==============================================================================
df <- data.frame(
  Depth = rep(c(10, 20, 40, 60, 80), each = 4),
  Tool = rep(c("DenovoFusion", "HMFtools", "GeneFuse", "FACTERA"), 5),
  
  # Data points
  Standard = c(3,3,2,0,  3,3,3,0,  3,3,3,0,  3,3,3,0,  3,3,3,0),
  LowVAF   = c(0,0,0,0,  1,0,0,0,  3,1,1,0,  3,2,1,0,  2,2,1,0),
  HighSNP  = c(2,0,0,0,  2,2,0,0,  2,2,0,0,  1,2,0,0,  2,2,0,0)
)

# Calculate Total Sensitivity (%)
df$Total <- round((df$Standard + df$LowVAF + df$HighSNP) / 9 * 100, 1)

# Set Factor Levels (Ensure DenovoFusion is at the top)
df$Tool <- factor(df$Tool, levels = rev(c("DenovoFusion", "HMFtools", "GeneFuse", "FACTERA")))
df$Depth <- factor(df$Depth, 
                   levels = c(10, 20, 40, 60, 80),      
                   labels = c("10x", "20x", "40x", "60x", "80x"))

# ==============================================================================
# 2. Plotting Function (Heatmap Style)
# ==============================================================================
plot_heatmap_svg <- function(data, y_var, title_text, is_percent=FALSE) {
  
  # Determine fill limits and label formatting
  if (is_percent) {
    fill_limits <- c(0, 100)
    label_format <- function(x) paste0(x, "%")
    fill_color_low <- "#FEE0D2"  # Very light red
    fill_color_high <- "#CB181D" # Deep red
  } else {
    fill_limits <- c(0, 3)
    label_format <- function(x) x
    # Blue scale for counts
    fill_color_low <- "#DEEBF7"  # Very light blue
    fill_color_high <- "#08519C" # Deep blue
  }
  
  p <- ggplot(data, aes(x = Depth, y = Tool, fill = .data[[y_var]])) +
    # The Heatmap Tiles
    geom_tile(color = "white", size = 0.5) + 
    
    # The Text Labels inside tiles
    geom_text(aes(label = label_format(.data[[y_var]])), 
              color = "black", size = 5, fontface = "bold") +
    
    # Color Scale
    scale_fill_gradient(low = fill_color_low, high = fill_color_high, 
                        limits = fill_limits, name = "") +
    
    # Labels (Title added here)
    labs(x = "Sequencing Depth", y = "", title = title_text) +
    
    # Theme: Minimal and Clean
    theme_minimal(base_size = 14) + 
    theme(
      legend.position = "none",
      axis.text = element_text(color="black", size=12),
      axis.text.y = element_text(face="bold", size=13),
      axis.title.x = element_text(face="bold", size=12, margin = margin(t = 8)),
      
      # [MODIFIED] Small Title Settings
      plot.title = element_text(hjust = 0.5, face="bold", size=12, margin = margin(b = 8)),
      
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

# ==============================================================================
# 3. Generate and Save SVGs (With Titles)
# ==============================================================================

# 1. Standard Group
p1 <- plot_heatmap_svg(df, "Standard", "Standard Group (AF=0.5)")
ggsave("FigureXA.svg", plot = p1, device = "svg", width = 6, height = 4)

# 2. Low-VAF Group
p2 <- plot_heatmap_svg(df, "LowVAF", "Low-VAF Group (AF=0.2)")
ggsave("FigureXB.svg", plot = p2, device = "svg", width = 6, height = 4)

# 3. High-SNP Group
p3 <- plot_heatmap_svg(df, "HighSNP", "High-SNP Group (5% Rate)")
ggsave("FigureXC.svg", plot = p3, device = "svg", width = 6, height = 4)

# 4. Total Sensitivity
p4 <- plot_heatmap_svg(df, "Total", "Total Sensitivity (Recall rate)", is_percent=TRUE)
ggsave("FigureXD.svg", plot = p4, device = "svg", width = 6, height = 4)

print("✅ All 4 Heatmap figures with small titles have been saved successfully.")