#---------------------------------------------------------------
# Set the working directory to the current folder
#---------------------------------------------------------------
setwd(getwd())

#---------------------------------------------------------------
# Load necessary libraries
#---------------------------------------------------------------
library(readxl)  # to read Excel files
library(dplyr)   # data manipulation
library(tidyr)   # data reshaping (pivot_wider)
library(tibble)  # convert between rownames and columns
library(MOFA2)
library(reticulate)
library(ggplot2)
library(pROC)

#---------------------------------------------------------------
# Load the five omics layer files
#---------------------------------------------------------------
clinical <- read_xlsx("TNBC_For.Integration_Clinical.xlsx")        # Clinical data: Patient_ID, Prognosis
expr     <- read_xlsx("TNBC_For.Integration_Transcriptomics.xlsx") # Transcriptomics: Gene expression (Log2CPM)
meth     <- read_xlsx("TNBC_For.Integration_Methylation.xlsx")     # DNA Methylation: Beta values
snv      <- read_xlsx("TNBC_For.Integration_SNV.xlsx")             # SNV: mutation status (0/1)
cnv      <- read_xlsx("TNBC_For.Integration_CNV.xlsx")             # CNV: copy number status (0/1)

#---------------------------------------------------------------
# Ensure Patient_ID is the first column in each dataset
#---------------------------------------------------------------
clinical <- clinical %>% select(Patient_ID, everything())
expr     <- expr %>% select(Patient_ID, everything())
meth     <- meth %>% select(Patient_ID, everything())
snv      <- snv %>% select(Patient_ID, everything())
cnv      <- cnv %>% select(Patient_ID, everything())

#---------------------------------------------------------------
# Convert each layer into a matrix with Patient_ID as rownames
# Clinical: only Prognosis column is retained
#---------------------------------------------------------------
clinical_mat <- clinical %>%
  select(Patient_ID, Prognosis) %>%  # keep only relevant columns
  column_to_rownames("Patient_ID")   # convert Patient_ID to rownames

# Transcriptomics: pivot from long to wide format (patient × gene)
expr_mat <- expr %>%
  pivot_wider(
    id_cols = Patient_ID,  # rows = patients
    names_from = Gene,     # columns = genes
    values_from = Log2CPM  # values = expression
  ) %>%
  column_to_rownames("Patient_ID") %>%
  as.matrix()  # convert to numeric matrix

# Methylation: pivot from long to wide format (patient × gene)
meth_mat <- meth %>%
  pivot_wider(
    id_cols = Patient_ID,
    names_from = Gene,
    values_from = Beta
  ) %>%
  column_to_rownames("Patient_ID") %>%
  as.matrix()

# SNV: pivot from long to wide format (patient × gene)
snv_mat <- snv %>%
  pivot_wider(
    id_cols = Patient_ID,
    names_from = Gene,
    values_from = Mutation,
    values_fn = max,        # collapse duplicates (0/1)
    values_fill = 0         # missing = not mutated
  ) %>%
  column_to_rownames("Patient_ID") %>%
  as.matrix()

# CNV: pivot from long to wide format (patient × gene)
cnv_mat <- cnv %>%
  pivot_wider(
    id_cols = Patient_ID,
    names_from = Gene,
    values_from = CNV
  ) %>%
  column_to_rownames("Patient_ID") %>%
  as.matrix()

#---------------------------------------------------------------
# Standardize rownames (patients) and ensure all matrices have the same patient order
#---------------------------------------------------------------
common_patients <- rownames(clinical_mat)  # take patient order from clinical layer
expr_mat <- expr_mat[common_patients, , drop = FALSE]
meth_mat <- meth_mat[common_patients, , drop = FALSE]
snv_mat  <- snv_mat[common_patients, , drop = FALSE]
cnv_mat  <- cnv_mat[common_patients, , drop = FALSE]

# At this point:
# - Each matrix has rows = patients, columns = features (genes or variables)
# - All matrices are aligned by patient, ready for integration

#---------------------------------------------------------------
# Create a list to hold four omics layers
#---------------------------------------------------------------
omics_list <- list(
  Transcriptome = expr_mat,
  Methylation   = meth_mat,
  SNV           = snv_mat,
  CNV           = cnv_mat
)

#---------------------------------------------------------------
# Transpose omics matrices for MOFA
# MOFA expects: features (genes/probes) x samples (patients)
#---------------------------------------------------------------
omics_list_t <- lapply(omics_list, function(x) {
  t(as.matrix(x))
})

#---------------------------------------------------------------
# Convert Clinical Prognosis to numeric (Bad = 1, Good = 0)
#---------------------------------------------------------------
prognosis_numeric <- ifelse(clinical_mat$Prognosis == "Bad", 1, 0)
names(prognosis_numeric) <- rownames(clinical_mat)

#---------------------------------------------------------------
# Create metadata (rows = samples)
#---------------------------------------------------------------
meta_df <- data.frame(
  Patient_ID = rownames(clinical_mat),
  Prognosis  = prognosis_numeric,
  row.names  = rownames(clinical_mat)
)

#---------------------------------------------------------------
# Create the MOFA object
#---------------------------------------------------------------
mofa_object <- create_mofa(
  omics_list_t,
  meta_data = meta_df
)

#---------------------------------------------------------------
# Get default MOFA options
#---------------------------------------------------------------
data_opts  <- get_default_data_options(mofa_object)
model_opts <- get_default_model_options(mofa_object)
train_opts <- get_default_training_options(mofa_object)

#---------------------------------------------------------------
# Set MOFA model options
#---------------------------------------------------------------
# Set number of latent factors (adjust based on data)
model_opts$num_factors <- 3

# View-specific likelihoods (defaults are usually fine)
model_opts$likelihoods <- c(
  Transcriptome = "gaussian",
  Methylation   = "gaussian",
  SNV           = "bernoulli",
  CNV           = "bernoulli"
)

#---------------------------------------------------------------
# Set MOFA training options
#---------------------------------------------------------------
# Make training explicit and verbose
train_opts$verbose <- TRUE

# Adjust convergence and iterations
train_opts$convergence_mode <- "medium"  # options: slow, medium, fast
train_opts$maxiter <- 1000               # maximum iterations

# Optional: verbosity for monitoring progress
train_opts$seed <- 42                    # reproducibility

#---------------------------------------------------------------
# Prepare the MOFA object for training
#---------------------------------------------------------------
mofa_object <- prepare_mofa(
  object           = mofa_object,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

#---------------------------------------------------------------
# Train the MOFA model
#---------------------------------------------------------------
# Adjust this path to where Python is installed
reticulate::use_python("C:/Users/zapan/AppData/Local/Programs/Python/Python310/python.exe", required = TRUE)
reticulate::py_run_string("import mofapy2")

mofa_object <- run_mofa(mofa_object)

#---------------------------------------------------------------
# Examine variance explained by each factor
#---------------------------------------------------------------
# Total variance explained per view
plot_variance_explained(mofa_object, plot_total = TRUE)

# Get values to plot in Prism
var_exp <- get_variance_explained(mofa_object)
r2_total <- var_exp$r2_total$group1
r2_total

# Variance explained per factor per view
plot_variance_explained(mofa_object, plot_total = FALSE)

#---------------------------------------------------------------
# Plot absolute loadings of the top 5 features of selected factors
#---------------------------------------------------------------

# Select layer and factor
layer <- "Methylation"
factor <- "Factor3"

# Get loadings
mofa_weights <- get_weights(mofa_object)
w <- mofa_weights[[layer]][, factor]

# Top 5 features by absolute value
top5 <- sort(abs(w), decreasing = TRUE)[1:5]

# Keep original signed values
top5_signed <- w[names(top5)]

# Absolute values
abs_vals <- abs(top5_signed)

# Proportional rescaling: max = 1
top5_rescaled <- abs_vals / max(abs_vals)

# Create data frame for plotting
df_plot <- data.frame(
  Feature = names(top5_rescaled),
  AbsLoading = abs(top5_rescaled)
)

# Keep the order as in top5_rescaled
df_plot$Feature <- factor(df_plot$Feature, levels = rev(names(top5_rescaled)))
df_plot$Highlight <- ifelse(df_plot$AbsLoading == max(df_plot$AbsLoading), "Top", "Other")

# Lollipop Plot 
ggplot(df_plot, aes(x = AbsLoading, y = Feature)) +
  geom_segment(aes(x = 0, xend = AbsLoading, y = Feature, yend = Feature, color = Highlight), size = 3) +
  geom_point(aes(color = Highlight), size = 4) +
  scale_color_manual(values = c(Top = "red", Other = "black")) +
  labs(x = "Absolute loading on Factor 3", y = NULL) +
  theme_classic(base_size = 30) +
  theme(legend.position = "none",
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black")
  )

#---------------------------------------------------------------
# Is the factor value significantly associated with Good vs Bad prognosis?
#---------------------------------------------------------------

# Extract factor values
factors <- get_factors(mofa_object, factors = "all")$group1

# Convert to data frame
factors_df <- as.data.frame(factors)

# Add clinical outcome
factors_df$Prognosis <- prognosis_numeric

head(factors_df)

### Non-parametric group comparison ###

wilcox_results <- lapply(colnames(factors_df)[grepl("Factor", colnames(factors_df))], function(f) {
  
  test <- wilcox.test(
    factors_df[[f]] ~ factors_df$Prognosis
  )
  
  data.frame(
    Factor = f,
    P_value = test$p.value
  )
})

wilcox_df <- do.call(rbind, wilcox_results)
wilcox_df$FDR <- p.adjust(wilcox_df$P_value, method = "BH")

wilcox_df

# Plot Factor 1
ggplot(factors_df, aes(x = as.factor(Prognosis), y = Factor1)) +
  geom_boxplot(fill = "grey90", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3) +
  scale_x_discrete(labels = c("Good", "Bad")) +
  labs(
    x = "Prognosis",
    y = "MOFA Factor 1"
  ) +
  annotate("text", 
           x = 1.5, 
           y = max(factors_df$Factor1) * 1.10,
           label = "Wilcoxon FDR = 5.0e-02",
           size = 7) +
  theme_classic(base_size = 26)

# Plot Factor 2
ggplot(factors_df, aes(x = as.factor(Prognosis), y = Factor2)) +
  geom_boxplot(fill = "grey90", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3) +
  scale_x_discrete(labels = c("Good", "Bad")) +
  labs(
    x = "Prognosis",
    y = "MOFA Factor 2"
  ) +
  annotate("text", 
           x = 1.5, 
           y = max(factors_df$Factor2) * 1.2,
           label = "Wilcoxon FDR = 1.2e-03",
           size = 7) +
  theme_classic(base_size = 26)

# Plot Factor 3
ggplot(factors_df, aes(x = as.factor(Prognosis), y = Factor3)) +
  geom_boxplot(fill = "grey90", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3) +
  scale_x_discrete(labels = c("Good", "Bad")) +
  labs(
    x = "Prognosis",
    y = "MOFA Factor 3"
  ) +
  annotate("text", 
           x = 1.5, 
           y = max(factors_df$Factor3) * 1.3,
           label = "Wilcoxon FDR = 3.4e-06",
           size = 7) +
  theme_classic(base_size = 26)

### ROC (predictive performance) ###

# Factor 1
roc_f1 <- roc(
  response = factors_df$Prognosis,
  predictor = factors_df$Factor1
)

# Factor 2
roc_f2 <- roc(
  response = factors_df$Prognosis,
  predictor = factors_df$Factor2
)

# Factor 3
roc_f3 <- roc(
  response = factors_df$Prognosis,
  predictor = factors_df$Factor3
)

# Get ROC values
auc(roc_f1)
auc(roc_f2)
auc(roc_f3)

# Plot ROC results
plot(roc_f1, 
     col = "grey50", 
     lwd = 2, 
     main = "ROC curves for MOFA factors",
     cex.axis = 1.5,  # increases x/y tick labels
     cex.lab = 1.5    # increases axis titles
)

lines(roc_f2, col = "blue", lwd = 2)
lines(roc_f3, col = "red", lwd = 2)

legend("bottomright", 
       legend = c("Factor1 (AUC=0.66)", "Factor2 (AUC=0.77)", "Factor3 (AUC=0.88)"),
       col = c("grey50", "blue", "red"), lwd = 2,
       cex = 1.3  # increases legend text size
)

# AUC ~ 0.5 → no predictive value
# AUC > 0.7 → good
# AUC > 0.8 → strong

#---------------------------------------------------------------
# Functional Enrichment
#---------------------------------------------------------------
# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Function run_enrichment
run_enrichment <- function(mofa_object, layer, factor, top_n = 300) {
  
  # Extract weights
  w <- get_weights(mofa_object)[[layer]][, factor]
  
  # Rank by absolute contribution
  top_genes <- names(sort(abs(w), decreasing = TRUE))[1:top_n]
  
  # Convert gene symbols to Entrez IDs
  gene_map <- bitr(
    top_genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  # GO Biological Process enrichment
  ego <- enrichGO(
    gene          = gene_map$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(ego)
}

### Factor 2 enrichment ###

ego_F2 <- run_enrichment(
  mofa_object = mofa_object,
  layer = "Transcriptome",
  factor = "Factor2",
  top_n = 300
)

# Plot
dotplot(ego_F2, showCategory = 5) +
  ggtitle("Transcriptome Factor 2") +
  theme_classic(base_size = 16)

### Factor 3 enrichment ###

ego_F3 <- run_enrichment(
  mofa_object = mofa_object,
  layer = "Methylation",
  factor = "Factor3",
  top_n = 300
)

# Plot
dotplot(ego_F3, showCategory = 5) +
  ggtitle("Methylation Factor 3") +
  theme_classic(base_size = 16)















































