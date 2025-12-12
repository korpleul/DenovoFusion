# This script created the function domains on gene EWSR1 and FLI1, and the relationship with breakpoints which divide the gene in
# different exon area, see Figure6B

library(ggplot2)
library(dplyr)

# -------------------------------
# 1. Domain structure data

# EWSR1 domains mapped to exons
ews_domains <- data.frame(
  domain = c("FET family domain", "EWS-specific RRM", "RNA recognition motif (RRM)",
             "RNA-binding superfamily", "Zinc finger domain"),
  start = c(1, 10, 12, 13, 16),   # exon start
  end   = c(9, 11, 12, 15, 17),   # exon end
  category = c("FET family domain", "EWS-specific RRM", "RNA recognition motif (RRM)",
               "RNA-binding superfamily", "Zinc finger domain")
)

# FLI1 domains mapped to exons
fli1_domains <- data.frame(
  domain = c("SAM/Pointed domain", "ETS DNA-binding domain"),
  start = c(1, 6),   # exon start
  end   = c(1, 11),  # exon end
  category = c("SAM/Pointed domain", "ETS DNA-binding domain")
)

# -------------------------------
# 2. Fusion sample information (breakpoint location)
samples <- data.frame(
  Sample = c("IC009", "IC015", "IC049", "IC086", "IC114"),
  ews_exon_break = c(11, 7.5, 8.5, 10.5, 7.5),     
  fli_exon_break = c(6.5, 5.5, 6.5, 6.5, 8.5)
)

# -------------------------------
# 3. Drawing functions
draw_fusion <- function(sample_id, ews_break, fli_break, y_pos) {
  scale_factor <- 35
  offset <- 400  # 调整 offset，使比例合理
  
  # Main gene lines
  ews_line <- geom_segment(aes(
    x = 0, 
    xend = max(ews_domains$end) * scale_factor + 5,  # buffer 5 units
    y = y_pos + 0.2, yend = y_pos + 0.2
  ), color = "#1f77b4", size = 2)
  
  fli_line <- geom_segment(aes(
    x = offset, 
    xend = offset + (max(fli1_domains$end) + 1) * scale_factor,  # buffer +1 exon
    y = y_pos - 0.2, yend = y_pos - 0.2
  ), color = "#1f77b4", size = 2)
  
  
  # Domain blocks
  domain_ews <- geom_rect(data = ews_domains,
                          aes(xmin = start * scale_factor,
                              xmax = end * scale_factor, 
                              ymin = y_pos + 0.12, ymax = y_pos + 0.28, fill = category),
                          inherit.aes = FALSE)
  
  domain_fli <- geom_rect(data = fli1_domains,
                          aes(xmin = offset + start * scale_factor,
                              xmax = offset + pmax(end, start + 0.5) * scale_factor,
                              ymin = y_pos - 0.28, ymax = y_pos - 0.12, fill = category),
                          inherit.aes = FALSE)
  
  
  # Breakpoints
  bp_ews <- geom_segment(aes(x = ews_break * scale_factor,
                             xend = ews_break * scale_factor,
                             y = y_pos + 0.12,
                             yend = y_pos + 0.28),
                         color = "red", linetype = "dashed", size = 2)
  bp_fli <- geom_segment(aes(x = offset + fli_break * scale_factor,
                             xend = offset + fli_break * scale_factor,
                             y = y_pos - 0.28,
                             yend = y_pos - 0.12),
                         color = "blue", linetype = "dashed", size = 2)
  
  # Breakpoint labels
  label_ews <- annotate("text", x = -50, y = y_pos + 0.2, 
                        label = "EWSR1", color = "red", size = 8, fontface = "bold")
  label_fli <- annotate("text", x = -50, y = y_pos - 0.1, 
                        label = "FLI1", color = "blue", size = 8, fontface = "bold")
  
  # Exon labels
  exon_label_ews <- annotate("text", x = ews_break * scale_factor, y = y_pos + 0.05,
                             label = paste0("exon ", floor(ews_break)), color = "red", size = 8, fontface = "bold")
  exon_label_fli <- annotate("text", x = offset + fli_break * scale_factor, y = y_pos - 0.03,
                             label = paste0("exon ", floor(fli_break)), color = "blue", size = 8, fontface = "bold")
  
  # Sample label
  label_sample <- annotate("text", x = -90, y = y_pos + 0.6, label = sample_id, hjust = 0, size = 8, fontface = "bold")
  
  list(ews_line, fli_line, domain_ews, domain_fli, bp_ews, bp_fli, 
       label_ews, label_fli, label_sample, exon_label_ews, exon_label_fli)
}

# -------------------------------
# 4. Main drawing
p <- ggplot() +
  theme_minimal(base_family = "Arial") +
  theme(text = element_text(face = "bold"))

for (i in 1:nrow(samples)) {
  p <- p + draw_fusion(
    sample_id = samples$Sample[i],
    ews_break = samples$ews_exon_break[i],
    fli_break = samples$fli_exon_break[i],
    y_pos = 1 - (i - 1) * 1.2
  )
}

p <- p +
  scale_fill_manual(
    name = "Category",
    values = c(
      "FET family domain" = "#2ca02c",           
      "EWS-specific RRM" = "#9467bd",             
      "RNA-binding superfamily" = "#8c564b",     
      "Zinc finger domain" = "#17becf",           
      "Nucleotide-binding domain" = "#7f7f7f",    
      "SAM/Pointed domain" = "#8c8c8c",          
      "ETS DNA-binding domain" = "#ff7f0e",      
      "Other" = "#F5F5F5"  
    )
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.position = "bottom",
    text = element_text(face = "bold")
  ) +
  ggtitle("EWSR1-FLI1 Fusion Domain Structure") +
  theme(plot.title = element_text(size = 24, hjust = 0.5, margin = margin(t = 10, b = 10)))

p
