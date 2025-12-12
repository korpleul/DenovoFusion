# This script create the frequency of different exon bindings in EWSR1-FLI1
# See Figure6A 

library(circlize)
library(readxl)
library(dplyr)
library(ggplot2)
library(grid)

# input the data from different sheets
df = read_excel("breakpoints.xlsx", sheet = "EWSR1-FLI1")

# Prepare the data
link_df <- df %>%
  count(exon_up1, exon_down2) %>%
  rename(from = exon_up1, to = exon_down2, value = n)

# Convert to character type (needed by circlize)
link_df$from <- paste0("EWSR1_exon", link_df$from)
link_df$to <- paste0("FLI1_exon", link_df$to)

# Define colors
sectors <- unique(c(link_df$from, link_df$to))
grid.col <- c(
  "EWSR1_exon7" = "#1f77b4",  
  "EWSR1_exon8" = "#ff7f0e",  
  "EWSR1_exon10" = "#2ca02c", 
  "EWSR1_exon11" = "#d62728", 
  "FLI1_exon5" = "#9467bd",   
  "FLI1_exon6" = "#8c564b",   
  "FLI1_exon8" = "#e377c2"    
)

# Output the figure
svg("Exon_intron_EWSR1-FLI1_chord.svg", width = 10, height = 10)

chordDiagram(link_df, grid.col = grid.col, transparency = 0.25, annotationTrack = "grid", preAllocateTracks = 1)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim), y = ylim[1] + 0.2,
      labels = gsub(".*_exon", "exon", sector.name),
      facing = "bending",     
      niceFacing = TRUE,
      adj = c(0.5, 0.5),
      cex = 1.6
    )
  },
  bg.border = NA
)

center_x <- 0.5
center_y <- 0.5
radius <- 0.5 

grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(center_y, center_y), "npc"),
           gp = gpar(col = "black", lwd = 1, lty = "dashed"))

title("Distribution of Exon-Exon Pairings in EWSR1-FLI1", line = -2, cex.main = 2.4)

dev.off()
