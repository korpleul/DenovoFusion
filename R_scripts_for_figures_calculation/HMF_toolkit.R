# 2. scatter chart
# Figure3C comparison between simple and complex SVs
set.seed(123)
simple_var <- c(32,46,25,33,44,38,32,52,80,31,
                56,35,43,86,33,41,55,46,68,95,
                37,106,85,52,41,42,43,42,80,37,
                62,73,60,33,37,30,65,40,42,39,
                43,35,40,96,49,109,41,32,38,48,
                31,49,35,27,25,56,45,51,61,63,
                39,36,55,56,42,51,39,38,38,38,
                27,38,39,36,87,41,72,68,40,47,
                34,60,38,74,34,54,45,62,46,30,
                81,36,41,44,43,52,43,42,41,50)  # Simple variations

complex_var <- c(108,86,101,104,97,85,89,98,109,106,
                 89,94,71,81,95,93,100,120,86,100,
                 100,102,79,102,85,100,97,91,113,93,
                 90,84,111,94,101,101,100,103,95,72,
                 89,92,100,104,91,95,100,81,98,92,
                 91,96,80,90,91,75,93,106,101,79,
                 88,97,92,101,82,116,80,73,107,87,
                 100,80,95,81,99,90,98,75,84,103,
                 86,95,73,82,92,103,73,99,91,81,
                 91,74,94,92,93,94,92,106,92,91)  # Complex variations

# Create the data frame
data <- data.frame(Simple = simple_var, Complex = complex_var)

# Load ggplot2
library(ggplot2)

# Create the scatter plot with a single color for points
ggplot(data, aes(x = Simple, y = Complex)) +
  geom_point(color = "#1F77B4", size = 1.5) +  # Single color for all points
  labs(x = "Simple Variations per Sample", y = "Complex Variations per Sample") +
  scale_x_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 40)) +  # Set x-axis limits and breaks
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 40)) +  # Set y-axis limits and breaks
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.title.x = element_text(size = 25),  # Change x-axis title size
        axis.title.y = element_text(size = 25),  # Change y-axis title size
        plot.title = element_text(size = 16),  # Change plot title size
        axis.text.x = element_text(size = 25),  # Change x-axis text size
        axis.text.y = element_text(size = 25)   # Change y-axis text size
  )  # Adding arrows to axes



# 3. Histogram
# Figure3D deviation in first 20 samples
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Sample names
samples <- c("IC009T", "IC015T", "IC024T", "IC034T", "IC044T",
             "IC046T", "IC049T", "IC053T", "IC054T", "IC057T",
             "IC058T", "IC066T", "IC067T", "IC071T", "IC076T",
             "IC077T", "IC080T", "IC082T", "IC086T", "IC092T")

# Replace these values with your actual percentages
data <- data.frame(
  Sample = rep(samples, each = 6),
  Type = rep(c("Deletion", "Incomplete", "Complex", "LINE", "2-break reciprocal", "2-break other"), times = 20),
  Count = c(1,1,108,13,2,16,
            18,1,86,8,2,18,
            3,2,101,8,1,11,
            3,1,104,10,1,19,
            15,2,97,7,0,21,
            4,1,85,7,3,22,
            1,1,89,7,2,22,
            23,1,98,9,1,19,
            48,2,109,10,3,15,
            11,1,106,3,2,15,
            34,2,89,7,2,13,
            8,2,94,6,2,18,
            18,1,71,4,2,18,
            36,1,81,8,1,18,
            3,0,95,6,3,21,
            4,0,93,8,0,25,
            22,2,100,11,1,21,
            2,1,120,10,6,28,
            36,1,86,13,3,16,
            57,3,100,10,3,22)
)

# Normalize percentages so that they sum to 100 for each sample
data <- data %>%
  group_by(Sample) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Define colors for each Type
colors <- c(
  "Deletion" = "#D55E00",           # 深橙红，NEJM中常用强调色（代替鲜红）
  "Incomplete" = "#0072B2",         # 深海军蓝，稳重专业（代替亮蓝）
  "Complex" = "#E69F00",            # 深金橙色（代替亮橙）
  "LINE" = "#009E73",               # 稳重深绿（代替浅绿）
  "2-break reciprocal" = "#CC79A7", # 柔和紫红（代替亮紫）
  "2-break other" = "#F0E442"       # 柠檬黄，明亮但不刺眼
)

# Reverse the Sample factor levels
data$Sample <- factor(data$Sample, levels = rev(unique(data$Sample)))

# Proceed with the plot
p <- ggplot(data, aes(x = Sample, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "fill") +  # Use position = "fill" to stack percentages
  labs(x = "Samples", y = "Percentage of SVs") +  # Set y label
  theme_classic() +
  coord_flip() +  # Flip the coordinates to rotate the chart
  scale_y_reverse() +  # Reverse the y-axis to start from top
  scale_fill_manual(values = colors) +  # Set custom colors
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),  # Adjust x-axis labels for better visibility
        axis.ticks.x = element_blank(),  # Remove x-axis ticks 
        axis.title.y = element_blank(),  # Remove y-axis title
        panel.grid.major.y = element_blank(),  # Remove major grid lines for y-axis
        panel.grid.minor.y = element_blank(),  # Remove minor grid lines for y-axis
        axis.line.y = element_blank(),  # Remove y-axis line
        axis.line.x = element_blank(),  # Remove x-axis line
        legend.position = "none",  # Remove the legend from the plot
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +  # Adjust margins if needed
  scale_y_continuous(labels = scales::percent)  # Format y-axis as percentage

# Save the plot with specified width and height
ggsave(filename = "HMF_SV.svg", plot = p, width = 3, height = 6)


# build a figure legend
# Load required libraries
library(ggplot2)
library(grid)

# Define the color palette
colors <- c(
  "Deletion" = "#D55E00",           # 深橙红，NEJM中常用强调色（代替鲜红）
  "Incomplete" = "#0072B2",         # 深海军蓝，稳重专业（代替亮蓝）
  "Complex" = "#E69F00",            # 深金橙色（代替亮橙）
  "LINE" = "#009E73",               # 稳重深绿（代替浅绿）
  "2-break reciprocal" = "#CC79A7", # 柔和紫红（代替亮紫）
  "2-break other" = "#F0E442"       # 柠檬黄，明亮但不刺眼
)

# Create a dummy data frame just to generate the legend
legend_data <- data.frame(
  Type = c("Deletion", "Incomplete", "Complex", "LINE", 
           "2-break reciprocal", "2-break other"),
  Count = c(1, 1, 1, 1, 1, 1)
)

# Generate the separate legend plot
legend_plot <- ggplot(legend_data, aes(x = 1, y = Count, fill = Type)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = colors) + 
  theme_void() +  # Remove axes and labels for a clean legend
  theme(legend.position = "right",  # Position the legend on the right
        legend.title = element_blank(),  # Remove the legend title
        legend.text = element_text(size = 12))  # Set text size for better readability

# Extract only the legend using cowplot
library(cowplot)
legend <- get_legend(legend_plot)

# 新建一个绘图设备（例如PNG），设置背景透明
png("HMF_figure_legend.png", width = 400, height = 200, bg = "transparent", res = 150)

# Display the legend
grid.newpage()
grid.draw(legend)

dev.off()
