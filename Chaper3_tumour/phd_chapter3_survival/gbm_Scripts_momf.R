library(MOMF)

data_dir <- '/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_CGGA/'


### Bulk data prep
raw_counts_325 <-  read.delim(paste0(data_dir,"CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt"),check.names = FALSE)
raw_counts_693 <-  read.delim(paste0(data_dir,"CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt"),check.names = FALSE)

raw_counts_325[is.na(raw_counts_325)] <- 0
raw_counts_693[is.na(raw_counts_693)] <- 0
raw_counts_325 <- raw_counts_325[raw_counts_325$gene_name!='',]
raw_counts_325 <- raw_counts_325[!duplicated(raw_counts_325$gene_name),]
raw_counts_693 <- raw_counts_693[raw_counts_693$gene_name!='',]
raw_counts_693 <- raw_counts_693[!duplicated(raw_counts_693$gene_name),]

rownames(raw_counts_325) <- raw_counts_325$gene_name
rownames(raw_counts_693) <- raw_counts_693$gene_name


raw_counts_693_filt <- as.data.frame(t(raw_counts_693[-1]))
raw_counts_325_filt <- as.data.frame(t(raw_counts_325[-1]))
rownames(comb_clin_filt) <- comb_clin_filt$CGGA_ID


raw_counts_693_filt <- raw_counts_693_filt[str_sub(rownames(raw_counts_693_filt)) %in% row.names(comb_clin_filt), ]
raw_counts_325_filt <- raw_counts_325_filt[str_sub(row.names(raw_counts_325_filt)) %in% row.names(comb_clin_filt), ]


bulk_counts <- rbind(raw_counts_693_filt,raw_counts_325_filt)
bulk_counts <-t(bulk_counts)
### sc data prep 
suppressMessages(library(rhdf5))

DATA.PATH <- file.path('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_survival', "data", "processed_data.h5")

# Read matrix
barcodes <- h5read(DATA.PATH, "expression_matrix/barcodes")
features <- h5read(DATA.PATH, "expression_matrix/features")
data <- h5read(DATA.PATH, "expression_matrix/data")
indices <- h5read(DATA.PATH, "expression_matrix/indices")
indptr <- h5read(DATA.PATH, "expression_matrix/indptr")
length(features)

#expr.mtx <- Matrix::sparseMatrix(
#  i=indices, p=indptr, x=as.numeric(data),
#  dimnames=list(as.character(features), as.character(barcodes)),
#  index1=FALSE, repr="C")

# For CSR format from Python, transpose to get CSC in R
expr.mtx <- Matrix::sparseMatrix(
  i = indices, 
  p = indptr, 
  x = as.numeric(data),
  dims = c(length(features), length(barcodes)),
  dimnames = list(as.character(features), as.character(barcodes)),
  index1 = FALSE, 
  repr = "C"
)

common_genes <- intersect(rownames(expr.mtx), rownames(bulk_counts))


expr_dense<- expr.mtx[common_genes,]
bulk_counts<- bulk_counts[common_genes,]

###
expr_dense <- as.matrix(expr_dense)

### Get cell annotations
sc_cell_type<- read.delim('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_survival/cell_type.csv',
                          sep=',',row.names = 'CellID')
sc_anno <- as.character(sc_cell_type$custom_annotation)

### compute the cell type specific expression level as reference
priorU <- momf.computeRef(expr_dense, sc_anno)

### create the gene list for MOMF 
GList <- list(X1 = t(expr_dense), X2 = t(bulk_counts))

### run MOMF
momf_res <- momf.fit(DataX = GList, DataPriorU = priorU, method = "KL", rho = 2, num_iter = 100)

### output the cell type proportions
cell_prop <- momf_res$cell.prop
### Save 
write.csv(cell_prop,'/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_survival/gbm_OUTPUT_momf_cell_prop.csv', 
          row.names = TRUE)

hist(cell_prop$oligo)

heatmap(cell_prop)

df <- data.frame(
  Sample = rownames(cell_prop),
  Oligodendrocytes = cell_prop[, "Oligodendrocytes"]
)
df <- df[order(df$Oligodendrocytes, decreasing = TRUE), ]


Sol <- ggplot(df, aes(x = Oligodendrocytes)) +
  geom_histogram(
    bins = 30,
    fill = "#4F81BD",
    color = "white",
    alpha = 0.8
  ) +
  theme_minimal(base_size = 8) +
  labs(x = "Oligodendrocyte proportion", y = "Frequency")+theme()

ggsave(
  filename = paste0(fig_outdir, "SM_ol_dist.png"),
  plot = Sol,
  device = "png",
  width = 2.5,
  height = 2,
  units = "in",
  dpi = 300
)

Sol <- ggplot(comb_clin_filt, 
              aes(x = Sd._risk_score, y = ol_proportions, fill = Sd._risk_score)) +
  geom_boxplot(alpha = 0.7) +
  labs(y = "OLs proportions", x = '', fill = "Sd. risk score") +
  scale_fill_manual(values = c("High" = "#D55E00", "Low" = "#0072B2")) +
  theme_minimal(base_size = 8) +
  theme(legend.position = "top",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

Sol<- ggplot(comb_clin_filt, aes(x = Sd._risk_score, y = ol_proportions, fill = Sd._risk_score)) +
  geom_boxplot(alpha = 0.7) +labs(y = "OLs proportions", x='',fill="Sd. risk score")+
  theme_minimal(base_size = 8)+theme()

ggsave(
  filename = paste0(fig_outdir, "SM_ol_grouping.png"),
  plot = Sol,
  device = "png",
  width = 2.5,
  height = 1.5,
  units = "in",
  dpi = 300
)

############

# Test if OL proportion correlates with risk score
ol_proportions <- cell_prop[, 3]
names(ol_proportions) <- rownames(cell_prop)

comb_clin_filt$ol_proportions <- ol_proportions[rownames(comb_clin_filt)]

# Correlation test
cor_test <- cor.test(ol_proportions, comb_clin_filt$z_score, method = "spearman")

# Cox model with OL proportion
cox_ol_prop <- coxph(surv_obj ~ ol_proportions)

summary(cox_ol_prop)

# Cox model: risk score + OL proportion
cox_combined <- coxph(surv_obj ~ z_score + ol_proportions,data= comb_clin_filt)
summary(cox_combined)

median_ol <- median(ol_proportions[ol_proportions > 0])
comb_clin_filt$ol_group <- ifelse(ol_proportions > median_ol, "High OL", "Low OL")

fit_ol <- survfit(surv_obj ~ ol_group , data = comb_clin_filt)


sol_kp<-ggsurvplot(
  fit_ol,
  data = comb_clin_filt,pval = TRUE,  pval.size = 3, pval.coord = c(x = 2000, y = 0.3),  
  palette = c("#D55E00", "#0072B2"),
  xlab = "Days",
  ylab = "Overall survival",
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  break.x.by = 500,
  font.main = c(8, "plain"),       # size, style
  font.x = c(8, "plain"),
  font.y = c(8, "plain"),
  font.tickslab = c(8, "plain"),
  font.legend = c(8, "plain"),
  risk.table.fontsize = 8 / 2.845,   # roughly 8 pt
  tables.theme = theme(
    plot.title = element_text(size = 8),
    text = element_text(size = 8, family = "sans"),
    axis.text.x = element_text(size =8, family = "sans"),
    axis.text.y = element_text(size = 8, family = "sans"),
    axis.title.x=  element_text(size = 8, family = "sans"),
    axis.title.y=  element_text(size = 8, family = "sans"),
    axis.title= element_text(size = 8, family = "sans"),
    axis.text.x.top = element_text(size = 8, family = "sans")))

ggsave(
  filename = paste0(fig_outdir, "SM_ol_kp_title.png"),
  plot = sol_kp$plot,
  device = "png",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)


ggsave(
  filename = paste0(fig_outdir, "SM_ol_table.png"),
  plot = sol_kp$table,
  device = "png",
  width = 3,
  height = 1,
  units = "in",
  dpi = 300
)

# Does signature work in both groups?
survdiff(surv_obj ~ Sd._risk_score, subset = ol_group == "High OL",data= comb_clin_filt)
survdiff(surv_obj ~ Sd._risk_score, subset = ol_group == "Low OL",data= comb_clin_filt)

###
