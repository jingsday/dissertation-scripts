# Assuming you have your expression matrix and survival data
set.seed(123)  # For reproducibility

# Your original 58 genes

# All genes in your expression matrix
all_genes <- colnames(RNA_combined_filt)

# Remove your 58 genes from the pool (optional, more conservative)
available_genes <- setdiff(all_genes, coef_data$variable)

# Generate 1000 random signatures
n_permutations <- 1000
random_scores <- matrix(NA, nrow = nrow(RNA_combined_filt), 
                        ncol = n_permutations)

for (i in seq_len(n_permutations)) {
  # Randomly sample 58 genes
  random_genes <- sample(available_genes, size = length(ol_genes), replace = FALSE)
  
  # Subset expression matrix to random genes
  random_expr <- RNA_combined_filt[, random_genes]
  
  # Compute random risk score
  # (example: mean expression per sample; adapt if using dot product)
  random_score <- rowMeans(random_expr)
  
  # Standardize (z-score)
  random_scores[, i] <- scale(random_score)
}


library(survival)
library(survcomp)

comb_clin_filt <- merge(comb_clin_filt, 
                        comb_clin["z_score"], 
                        by = "row.names", 
                        all.x = TRUE)


# Your original signature's C-index
ol_cindex <- concordance.index(
  x = comb_clin_filt$z_score,
  surv.time = comb_clin_filt$OS,
  surv.event = comb_clin_filt$Censor..alive.0..dead.1. == "1"
)$c.index

# C-index for each random signature
random_cindices <- numeric(n_permutations)

for(i in 1:n_permutations) {
  random_cindices[i] <- concordance.index(
    x = random_scores[, i],
    surv.time = comb_clin_filt$OS,
    surv.event = comb_clin_filt$Censor..alive.0..dead.1. == "1"
  )$c.index
}


# Calculate empirical p-value
p_value <- mean(random_cindices >= ol_cindex)

# Summary statistics
cat("OL signature C-index:", round(ol_cindex, 3), "\n")
cat("Random signatures C-index (mean ± SD):", 
    round(mean(random_cindices), 3), "±", 
    round(sd(random_cindices), 3), "\n")
cat("Empirical p-value:", p_value, "\n")


library(ggplot2)

# Create dataframe
df <- data.frame(
  C_index = random_cindices,
  Type = "Random Signatures"
)

# Add your OL signature
ol_data <- data.frame(
  C_index = ol_cindex,
  Type = "OL Signature"
)

# Combine
plot_data <- rbind(df, ol_data)

# Create figure
ggplot() +
  # Histogram of random signatures
  geom_histogram(data = df, aes(x = C_index), 
                 bins = 50, fill = "gray70", color = "black", alpha = 0.7) +
  # Vertical line for OL signature
  geom_vline(xintercept = ol_cindex, color = "#D32F2F", 
             linewidth = 1.5, linetype = "solid") +
  # Add annotation
  annotate("text", x = ol_cindex - 0.05, y = 80,
           label = paste0("OL Signature\nC-index = 0.838"),
           color = "#D32F2F", fontface = "bold", hjust = 1, size = 5) +
  annotate("text", x = 0.52, y = 80,
           label = paste0("Random Signatures\nMean = 0.484 ± 0.028\nn = 1,000"),
           color = "gray30", fontface = "bold", size = 4.5) +
  # Styling
  labs(x = "Concordance Index (C-index)", 
       y = "Frequency",
       title = "OL Signature Outperforms Random Gene Sets") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  scale_x_continuous(limits = c(0.4, 0.9), breaks = seq(0.4, 0.9, 0.1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

library(ggplot2)

t<- ggplot() +
  # Histogram of random signatures
  geom_histogram(data = df, aes(x = C_index), 
                 bins = 50, fill = "gray70", color = "black", alpha = 0.7) +
  
  # Vertical line for OL signature
  geom_vline(xintercept = ol_cindex, color = "#D32F2F", 
             linewidth = 1.2, linetype = "solid") +
  
  # Add annotations
  annotate("text", x = ol_cindex - 0.03, y = 65,
           label = paste0("OL Signature\nC-index = ", sprintf("%.3f", ol_cindex)),
           color = "#D32F2F",  hjust = 1, size = 8/2.8) +  # ~8 pt
  
  annotate("text", x = 0.7, y = 130,
           label = paste0("Random Signatures\nMean = ",
                          sprintf("%.3f ± %.3f", mean(random_cindices), sd(random_cindices)),
                          "\nn = ", n_permutations),
           color = "black", size = 8/2.8) +  # ~8 pt
  
  # Labels
  labs(x = "Concordance Index (C-index)", 
       y = "Frequency"# title = "OL Signature Outperforms Random Gene Sets"
       ) +
  
  # Theme
  theme_classic(base_size = 8) +
  theme(#plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    axis.title = element_text( size = 8),
    axis.text  = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  
  # X-axis limits and breaks
  scale_x_continuous(limits = c(0.4, 0.9), breaks = seq(0.4, 0.9, 0.1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

t


ggsave(
  filename = paste0(fig_outdir, "F3_specificity.png"),
  plot = t,
  device = "png",
  width = 3,
  height = 3,
  units = "in",
  dpi = 300
)
###############
# Test 100 random signatures and compare their HRs
set.seed(123)
n_tests <- 100
random_hrs <- numeric(n_tests)

for(i in 1:n_tests) {
  # Random 58 genes
  random_genes <- sample(rownames(expression_matrix), 58)
  random_score <- colMeans(log2(expression_matrix[random_genes, ] + 1))
  random_score_std <- scale(random_score)
  
  # Fit Cox model
  cox_model <- coxph(Surv(time, event) ~ random_score_std)
  random_hrs[i] <- exp(coef(cox_model))
}

# Your OL signature HR = 8.45
# Compare
mean(random_hrs > 8.45)  # Probability random signature does better
