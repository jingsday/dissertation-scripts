library(ComplexHeatmap)
library(circlize)
library(grid) # required for unit() and gpar()

patient_order <- order(comb_clin_filt$z_score)

# Reorder data
cell_prop_ordered <- cell_prop[patient_order, ]
cell_prop_t <- t(cell_prop_ordered)

# Extract OL proportions
ol_proportions <- cell_prop[, "Oligodendrocytes"]

# Create annotation bars with adjusted heights

###
ha_top <- HeatmapAnnotation(
  `Risk Score` = comb_clin_filt$z_score[patient_order],
  `Survival` = ifelse(comb_clin_filt$Censor..alive.0..dead.1.[patient_order] == 0, 
                      "Alive", "Dead"),
  
  col = list(
    `Risk Score` = colorRamp2(
      c(-2.5, 0, 2.5), 
      c("#0072B2", "white", "#D55E00")
    ),
    `Survival` = c("Alive" = "gray80", "Dead" = "black")
  ),
  
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  gap = unit(2, "mm"),
  show_legend = c(TRUE, TRUE),
  border = TRUE,
  
  annotation_legend_param = list(
    `Risk Score` = list(
      title = "Risk Score (z)",
      at = seq(-2.5, 2.5, by = 1),
      labels = seq(-2.5, 2.5, by = 1),
      legend_height = unit(4, "cm"),
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8)
    ),
    `Survival` = list(
      title = "Survival",
      legend_height = unit(2, "cm"),
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8)
    )
  )
)

ha_top

# Adjust color scale to show more variation (optional adjustments commented out)
cell_prop_t_adj <- cell_prop_t 

# Create heatmap with NO row clustering to show all cell types clearly
ht <- Heatmap(
  cell_prop_t,
  name = "Cell\nProportion",
  
  # Adjusted color scale - more sensitive to small values
  col = colorRamp2(
    seq(0, 1, length.out = 9),
    RColorBrewer::brewer.pal(9, "Reds")
  ),
  top_annotation = ha_top,
  
  # DISABLE clustering to show all rows
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  
  # Larger row labels
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  row_names_max_width = unit(4, "cm"),
  
  show_column_names = FALSE,
  column_title = "Recurrent GBM Patients (n=88, ordered by risk score)",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  
  # Better spacing
  row_gap = unit(0.5, "mm"),
  
  heatmap_legend_param = list(
    title = "Cell\nProportion",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(5, "cm"),
    direction = "vertical",
    at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0")
  ),
  
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.8),
  
  # Increase heatmap height and width
  height = unit(2, "inch"),
  width = unit(4, "inch")
)

draw(ht) 
fig_outdir <- '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_fig/'
png(
  filename = paste0(fig_outdir, "SM_heatmap.png"),
  width = 6.5,
  height = 4,
  units = "in",
  res = 300
)
draw(ht)
dev.off()


###### simpler version
####Create color palettes
risk_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
risk_breaks <- seq(min(comb_clin_filt$z_score), max( comb_clin_filt$z_score), length.out = 101)
risk_col <- risk_colors[cut( comb_clin_filt$z_score[patient_order], breaks = risk_breaks)]

group_col <- ifelse( comb_clin_filt$Sd._risk_score[patient_order] == "High", "#E91E63", "#4CAF50")

ol_colors <- colorRampPalette(c("white", "#1565C0"))(100)
ol_breaks <- seq(0, max(ol_proportions), length.out = 101)
ol_col <- ol_colors[cut(ol_proportions[patient_order], breaks = ol_breaks)]

# Setup layout
layout(matrix(c(1,1,2,2,2,2,2,2), nrow = 8, ncol = 1))
par(mar = c(0, 8, 2, 2))

# Top annotations
plot.new()
plot.window(xlim = c(0, ncol(cell_prop)), ylim = c(0, 3))

# Risk score bar
rect(0:(ncol(cell_prop)-1), 2, 1:ncol(cell_prop), 2.3, 
     col = risk_col, border = NA)
text(-2, 2.15, "Risk Score", adj = 1, font = 2, xpd = TRUE)

# Risk group bar
rect(0:(ncol(cell_prop)-1), 1.6, 1:ncol(cell_prop), 1.9, 
     col = group_col, border = NA)
text(-2, 1.75, "Risk Group", adj = 1, font = 2, xpd = TRUE)

# OL proportion bar
rect(0:(ncol(cell_prop)-1), 1.2, 1:ncol(cell_prop), 1.5, 
     col = ol_col, border = NA)
text(-2, 1.35, "OL %", adj = 1, font = 2, xpd = TRUE)

# Main heatmap
par(mar = c(4, 8, 0, 2))
heatmap(
  cell_prop,
  Colv = NA,  # Don't cluster columns
  Rowv = NULL,  # Cluster rows
  scale = "none",
  col = colorRampPalette(c("white", "#FFF9C4", "#FFB74D", "#E64A19", "#B71C1C"))(50),
  labCol = rep("", ncol(cell_prop)),  # No column labels
  cexRow = 1.2,
  margins = c(5, 10),
  main = ""
)

