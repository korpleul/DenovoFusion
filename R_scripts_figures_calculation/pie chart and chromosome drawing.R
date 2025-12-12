
# pie chart for filter function
# Define the values for each category
total <- 46
no_valid_annotation <- 21
duplicates <- 3
long_gaps <- 4
small_fragments <- 4
homologs <- 1
min_support <- 7
remaining <- 6

# Create a vector for the categories
categories <- c("No Valid Annotation", "Duplicates", "Long Gaps", 
                "Small Fragments", "Homologs", "Min Support", "Remaining")

# Create a vector for the values (counts) corresponding to each category
counts <- c(no_valid_annotation, duplicates, long_gaps, small_fragments, 
            homologs, min_support, remaining)

# Calculate percentages
percentages <- round((counts / total) * 100, 1)

# Create labels with category names and percentages
labels <- paste(categories, "\n", percentages, "%", sep="")

svg("Figure4B.svg", width = 10, height = 10, bg = "transparent")
pie(counts, labels = labels, col = c(
  "#d62728",  
  "#1f77b4",  
  "#17becf",  
  "#ff7f0e",  
  "#7f7f7f",
  "#2ca02c",
  "#9467bd" 
),
    cex = 1.5)
dev.off()



# draw chromosome
library(Gviz)

# Define chromosomes and genome
chr22 <- "chr22"
chr11 <- "chr11"
chr_start <- 1
chr_end_22 <- 51304566  # Chromosome 22 length
chr_end_11 <- 134451414  # Chromosome 11 length

# Create the ideogram track for chromosome 22 (without the label and border)
chromosome_track_22 <- IdeogramTrack(genome = "hg38", chromosome = chr22, 
                                     border = NULL, 
                                     background.panel = "transparent", 
                                     background.title = "transparent",
                                     showId = FALSE)  # Remove label ("chr22")

# Create the ideogram track for chromosome 11 (without the label and border)
chromosome_track_11 <- IdeogramTrack(genome = "hg38", chromosome = chr11, 
                                     border = NULL, 
                                     background.panel = "transparent", 
                                     background.title = "transparent",
                                     showId = FALSE)  # Remove label ("chr11")

# Save the plot for chromosome 22 as SVG (vector format)
svg("chromosome_22_ideogram.svg", width = 10, height = 1, bg = "transparent")
plotTracks(chromosome_track_22,
           from = chr_start, to = chr_end_22,
           col.frame = "white",  # Ensure no border around the plot
           background.panel = "transparent",  # Transparent background
           background.title = "transparent",  # Transparent title background
           margin = 0)  # Remove margins
dev.off()

# Save the plot for chromosome 11 as SVG (vector format)
svg("chromosome_11_ideogram.svg", width = 10, height = 1,bg = "transparent")
plotTracks(chromosome_track_11,
           from = chr_start, to = chr_end_11,
           col.frame = "white",  # Ensure no border around the plot
           background.panel = "transparent",  # Transparent background
           background.title = "transparent",  # Transparent title background
           margin = 0)  # Remove margins
dev.off()

# draw a plot for breakpoints distribution
# 加载 ggplot2
library(ggplot2)

# 读取数据
data <- read.csv("/home/xinwei/Desktop/manuscript/breakpoints_information.csv")

# 
ews_range <- "EWSR1: 29,285,329 - 29,292,556 (chr22)"
fli_range <- "FLI1: 128,778,092 - 128,821,290 (chr11)"

svg("breakpoints.svg", width = 10, height = 6,bg = "transparent")
ggplot(data, aes(x = EWSR1, y = FLI1)) +
  geom_point(size = 2, alpha = 0.8) +  
  labs(title = "EWSR1-FLI1 breakpoints distribution", 
       x = "EWSR1 gene (chr22)", 
       y = "FLI1 gene置 (chr11)") +
  theme_classic() +  
  theme(axis.line = element_line(linewidth = 1, arrow = arrow(type = "closed")),  
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_continuous(limits = c(29268009, 29300525)) +  # X 轴范围
  scale_y_continuous(limits = c(128686535, 128813267)) +  # Y 轴范围
  annotate("text", x = 29268009 + 100, y = 128810000, label = ews_range, hjust = 0, size = 4, color = "black") +
  annotate("text", x = 29268009 + 100, y = 128805000, label = fli_range, hjust = 0, size = 4, color = "black") +
  annotate("rect", xmin = 29268009, xmax = 29270009, ymin = 128802000, ymax = 128813267, alpha = 0.2, fill = "grey")
dev.off()







