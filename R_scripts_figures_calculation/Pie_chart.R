# This file created the chart for fusion gene frequency from Curie group, HMFtools, 
# FACTERA, Genefuse and DenovoFusion, in manuscript Figure3B-E.
# 3d pie chart for fusion gene frequency from different methods
# Figure3B Frequency of EWSR1-FLI1 fusion in 100 samples
library(Cairo)
# Load necessary library
library(plotrix)

# Original from curie
# Define the data
counts <- c(74, 8, 1, 17)  # Frequencies of fusion and no fusion
labels <- c("EWSR1-FLI1", "EWSR1-ERG", "EWSR1-ETV1", "No Fusion")

# caculate the counts
total <- sum(counts)
labels_frac <- paste(labels, "\n", counts, "/", total, sep = "")

# Create a 3D pie chart
par(mar=c(0,0,0,0))  # Adjust margins if needed
NEJM_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#7f7f7f")
pie(counts, 
      labels = NA,  # Don't show default labels
      col = NEJM_colors,  # Custom colors
      radius = 1,  # Size of the pie
      border = NA)  # Adjust the viewing angle

text(-0.0,  0.6, labels_frac[1], cex = 2.5, font = 2, col = "black", pos = 1)
text(0.05, -0.85, labels_frac[2], cex = 2.5, font = 2, col = "black", pos = 1)
text(0.3, -0.5, labels_frac[3], cex = 2.5, font = 2, col = "black", pos = 1)
text( 0.6, -0.1, labels_frac[4], cex = 2.5, font = 2, col = "black", pos = 1)

# HMF toolkit
# Define the data
counts <- c(61, 4, 1, 34)  # Frequencies of fusion and no fusion
labels <- c("EWSR1-FLI1", "EWSR1-ERG", "EWSR1-ETV1", "No Fusion")

# caculate the counts
total <- sum(counts)
labels_frac <- paste(labels, "\n", counts, "/", total, sep = "")

# Create a 3D pie chart
par(mar=c(0,0,0,0))  # Adjust margins if needed
NEJM_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#7f7f7f")
pie(counts, 
      labels = NA,  # Don't show default labels
      col = NEJM_colors,  # Custom colors
      radius = 1,  # Size of the pie
      border = NA)  # Adjust the viewing angle

# Add custom percentage labels
pct <- round(counts / sum(counts) * 100, 1)  # Calculate percentages
text(-0.0, 0.7, labels_frac[1], cex = 2.5, font = 2, col = "black", pos = 1)
text(-0.4, -0.3, labels_frac[2], cex = 2.5, font = 2, col = "black", pos = 1)
text(-0.5, -0.7, labels_frac[3], cex = 2.5, font = 2, col = "black", pos = 1)
text(0.5, -0.2, labels_frac[4], cex = 2.5, font = 2, col = "black", pos = 1)

# FACTERA
# Define the data
counts <- c(0, 0, 0, 100)  # Frequencies of fusion and no fusion
labels <- c("EWSR1-FLI1", "EWSR1-ERG", "EWSR1-ETV1", "No Fusion")

# caculate the counts
total <- sum(counts)
labels_frac <- paste(labels, "\n", counts, "/", total, sep = "")

# Create a 3D pie chart
par(mar=c(0,0,0,0))  # Adjust margins if needed
NEJM_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#7f7f7f")
pie(counts, 
      labels = NA,  # Don't show default labels
      col = NEJM_colors,  # Custom colors
      radius = 1,  # Size of the pie
      border= NA ) 

# Add custom percentage labels
pct <- round(counts / sum(counts) * 100, 1)  # Calculate percentages
text(-0.0, 0.4, labels_frac[1], cex = 2.5, font = 2, col = "#7f7f7f", pos = 1)
text(-0.4, -0.3, labels_frac[2], cex = 2.5, font = 2, col = "#7f7f7f", pos = 1)
text(-0.4, -0.4, labels_frac[3], cex = 2.5, font = 2, col = "#7f7f7f", pos = 1)
text(0.0, 0.3, labels_frac[4], cex = 2.5, font = 2, col = "black", pos = 1)


# Genefuse
# Define the data
counts <- c(51, 4, 0, 45)  # Frequencies of fusion and no fusion
labels <- c("EWSR1-FLI1", "EWSR1-ERG", "EWSR1-ETV1", "No Fusion")

# caculate the counts
total <- sum(counts)
labels_frac <- paste(labels, "\n", counts, "/", total, sep = "")

# Create a 3D pie chart
par(mar=c(0,0,0,0))  # Adjust margins if needed
NEJM_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#7f7f7f")
pie(counts, 
    labels = NA,  # Don't show default labels
    col = NEJM_colors,  # Custom colors
    radius = 1,  # Size of the pie
    border= NA ) 

# Add custom percentage labels
pct <- round(counts / sum(counts) * 100, 1)  # Calculate percentages
text(-0.0, 0.6, labels_frac[1], cex = 2.5, font = 2, col = "black", pos = 1)
text(-0.6, -0.05, labels_frac[2], cex = 2.5, font = 2, col = "black", pos = 1)
text(0.4, -0.3, labels_frac[3], cex = 2.5, font = 2, col = "#7f7f7f", pos = 1)
text(0.4, -0.3, labels_frac[4], cex = 2.5, font = 2, col = "black", pos = 1)


# denovoFusion
# Define the data
counts <- c(27, 1, 0, 72)  # Frequencies of fusion and no fusion
labels <- c("EWSR1-FLI1", "EWSR1-ERG", "EWSR1-ETV1", "No Fusion")

# caculate the counts
total <- sum(counts)
labels_frac <- paste(labels, "\n", counts, "/", total, sep = "")

# Create a 3D pie chart
par(mar=c(0,0,0,0))  # Adjust margins if needed
NEJM_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#7f7f7f")
pie(counts, 
    labels = NA,  # Don't show default labels
    col = NEJM_colors,  # Custom colors
    radius = 1,  # Size of the pie
    border= NA ) 


# Add custom percentage labels
pct <- round(counts / sum(counts) * 100, 1)  # Calculate percentages
text(0.4, 0.4, labels_frac[1], cex = 2.5, font = 2, col = "black", pos = 1)
text(-0.1, 0.9, labels_frac[2], cex = 2.5, font = 2, col = "black", pos = 1)
text(0.3, -0.3, labels_frac[3], cex = 2.5, font = 2, col = "#7f7f7f", pos = 1)
text(-0.0, -0.2, labels_frac[4], cex = 2.5, font = 2, col = "black", pos = 1)

