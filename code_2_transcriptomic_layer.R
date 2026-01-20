# Set the working directory to the current directory
setwd(getwd())

# Load necessary libraries for data querying, manipulation, and analysis
library(TCGAbiolinks)
library(recount3)
library(SummarizedExperiment)
library(data.table)
library(sesame)
library(sesameData)
library(maftools)
library(limma)
library(edgeR)
library(writexl)
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(tidyr)


# Retrieve all available GDC projects
gdcprojects <- getGDCprojects()

# Get a summary of the TCGA-BRCA project
getProjectSummary('TCGA-BRCA')

#-----------------------------------------------
# Clinical Layer
#-----------------------------------------------

# Create a query to access clinical data from TCGA-BRCA project
query_clinical <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)

# Download the clinical data files based on the query above
# This step was done already in code_1
# GDCdownload(query_clinical)

# Prepare and load the downloaded clinical data into R as a SummarizedExperiment object
tcga_brca_clinical <- GDCprepare(query_clinical, summarizedExperiment = TRUE)

# List all clinical data tables available in the downloaded TCGA-BRCA clinical dataset
# There are 9 different clinical tables, each representing distinct clinical information
names(tcga_brca_clinical)

# Example: Retrieve and display the column names of the 'clinical_patient_brca' table
# This table contains key patient-level clinical variables
colnames(tcga_brca_clinical[[1]])

# Extract the patient clinical data
patient_data <- tcga_brca_clinical$clinical_patient_brca

# Check unique values for each receptor to understand how negatives are coded
unique(patient_data$er_status_by_ihc)
unique(patient_data$pr_status_by_ihc)
unique(patient_data$her2_status_by_ihc)

# Filter patients who are negative for ER, PR, and HER2 (triple negative)
tnbc_patients <- subset(patient_data,
                        er_status_by_ihc == "Negative" &
                          pr_status_by_ihc == "Negative" &
                          her2_status_by_ihc == "Negative"
)

# Count the number of triple-negative breast cancer patients
nrow(tnbc_patients)

# Get the patient barcodes or IDs
tnbc_patient_ids <- tnbc_patients$bcr_patient_barcode

# Output to excel (optional)
#tnbc_patient_ids_df <- as.data.frame(tnbc_patient_ids)
#write_xlsx(tnbc_patient_ids_df, "TNBC_Clinical_Patients.xlsx")

# Print how many triple negative patients were found
cat("Number of triple negative breast cancer patients:", length(tnbc_patient_ids), "\n")

# View patient IDs
print(tnbc_patient_ids)

# Result: 116 patients have TNBC.

#-----------------------------------------------
# Transcriptomic Layer
#-----------------------------------------------

# Create a query to access transcriptomic data from TCGA-BRCA project
query_transcriptomic <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# Download the transcriptomic data files based on the query above
GDCdownload(query_transcriptomic)

# Prepare and load the downloaded transcriptomic data into R as a SummarizedExperiment object
tcga_brca_transcriptomic <- GDCprepare(query_transcriptomic, summarizedExperiment = TRUE)

# Extract raw count matrix (genes x samples)
count_matrix <- assay(tcga_brca_transcriptomic)

# Extract sample barcodes
sample_barcodes <- colnames(count_matrix)

# Extract patient barcodes (first 12 chars)
sample_patients <- substr(sample_barcodes, 1, 12)

# Match TNBC patient barcodes
tnbc_sample_indices <- sample_patients %in% tnbc_patient_ids
count_matrix_tnbc <- count_matrix[, tnbc_sample_indices]

# Confirm dimensions
cat("Dimensions of TNBC count matrix: ", dim(count_matrix_tnbc), "\n")

# Extract barcodes of TNBC RNA-seq samples
tnbc_sample_barcodes <- colnames(count_matrix_tnbc)
tnbc_sample_patients <- substr(tnbc_sample_barcodes, 1, 12)

# Count RNA-seq samples per TNBC patient
table_tnbc_samples <- table(tnbc_sample_patients)

# Identify patients with duplicate RNA-seq samples
dup_patients <- names(table_tnbc_samples[table_tnbc_samples > 1])
cat("Duplicated patients:", dup_patients, "\n")

#-----------------------------------------------
# Resolve Duplicates by Comparing Quality
#-----------------------------------------------

# Store selected samples
selected_samples <- c()

# Loop through duplicated patients
for (pid in dup_patients) {
  
  # Get all samples for this patient
  sample_indices <- which(tnbc_sample_patients == pid)
  samples <- tnbc_sample_barcodes[sample_indices]
  counts_subset <- count_matrix_tnbc[, samples]
  
  # Compute metrics
  total_counts <- colSums(counts_subset)
  expressed_genes <- apply(counts_subset, 2, function(x) sum(x > 10))
  
  # Build metric table
  metrics_df <- data.frame(
    sample = samples,
    total_counts = total_counts,
    expressed_genes = expressed_genes,
    stringsAsFactors = FALSE
  )
  
  # Rank by expressed genes, then by total counts
  metrics_df <- metrics_df[order(-metrics_df$expressed_genes, -metrics_df$total_counts), ]
  
  # Keep the top sample
  selected_samples <- c(selected_samples, metrics_df$sample[1])
}

#-----------------------------------------------
# Add Unique Samples (No Duplicates)
#-----------------------------------------------

unique_samples <- tnbc_sample_barcodes[!(tnbc_sample_patients %in% dup_patients)]
final_selected_samples <- c(selected_samples, unique_samples)

# Final count matrix for TNBC
count_matrix_tnbc_filtered <- count_matrix_tnbc[, final_selected_samples]

# Extract unique sample barcodes
sample_barcodes_unique <- colnames(count_matrix_tnbc_filtered)

# Extract unique patient barcodes (first 12 chars)
sample_patients_ids <- substr(sample_barcodes_unique, 1, 12)

# Output to excel
tnbc_transcriptomic_ids_df <- as.data.frame(sample_patients_ids)
write_xlsx(tnbc_transcriptomic_ids_df, "TNBC_Transcriptomic_Patients.xlsx")

# Result: 115 TNBC patients have transcriptomic data.

#---------------------------------------------------------------
# Load the 60 TNBC patient list from Excel
#---------------------------------------------------------------
tnbc_final <- read_xlsx("TNBC_Patients_Final.xlsx", col_names = FALSE)

# Ensure patient barcodes are uppercase for consistency
final_patient_ids <- toupper(tnbc_final[[1]])  # Assuming first column contains patient IDs

#---------------------------------------------------------------
# Curated TNBC Patients
#---------------------------------------------------------------

# Select only your 60 curated TNBC patients
tnbc <- patient_data[patient_data$bcr_patient_barcode %in% final_patient_ids, ]

# Convert to a data.frame if it’s not already
tnbc_df <- as.data.frame(tnbc)

#---------------------------------------------------------------
# Prognosis definition
#---------------------------------------------------------------

# Ensure death_days_to is numeric
tnbc_df$death_days_to <- as.numeric(tnbc_df$death_days_to)

# Flag patients by each criterion

# 1. Early death (<5 years = 1825 days)
tnbc_df$bad_death <- ifelse(!is.na(tnbc_df$death_days_to) & tnbc_df$death_days_to < 1825, TRUE, FALSE)

# 2. Large/advanced tumor (T3 or T4, including subtypes like T4a)
tnbc_df$bad_tumor <- grepl("^T[34]", tnbc_df$ajcc_tumor_pathologic_pt, ignore.case = TRUE)

# 3. Extensive nodal involvement (N2 or N3, including subtypes like N3a)
tnbc_df$bad_nodes <- grepl("^N[23]", tnbc_df$ajcc_nodes_pathologic_pn, ignore.case = TRUE)

# 4. Advanced pathologic stage (Stage III or IV, including subtypes like IIIA)
tnbc_df$bad_stage <- grepl("Stage\\s*(III|IV)", tnbc_df$ajcc_pathologic_tumor_stage, ignore.case = TRUE)

# Combine tumor, nodes, and stage into one feature flag
tnbc_df$bad_prognosis <- tnbc_df$bad_death | tnbc_df$bad_tumor | tnbc_df$bad_nodes | tnbc_df$bad_stage

# Count how many patients have bad prognosis
num_bad <- sum(tnbc_df$bad_prognosis, na.rm = TRUE)
cat("Number of TNBC patients with bad prognosis:", num_bad, "\n")

# See patient IDs with bad prognosis
bad_patients <- tnbc_df$bcr_patient_barcode[tnbc_df$bad_prognosis]
print(bad_patients)

#---------------------------------------------------------------
# Differential Expression: Bad vs Good Prognosis TNBC Patients
#---------------------------------------------------------------

# Match transcriptomic samples to clinical data
prognosis_status <- tnbc_df$bad_prognosis[match(substr(colnames(count_matrix_tnbc_filtered), 1, 12), 
                                                tnbc_df$bcr_patient_barcode)]

# Filter samples with known prognosis
valid_indices <- !is.na(prognosis_status)
count_matrix_valid <- count_matrix_tnbc_filtered[, valid_indices]
prognosis_status <- prognosis_status[valid_indices]

# Convert logicals to factor labels
group_labels <- factor(ifelse(prognosis_status, "Bad", "Good"))

cat("Samples retained for analysis:", length(group_labels), "\n")
table(group_labels)

#---------------------------------------------------------------
# Create DGEList and Normalize
#---------------------------------------------------------------

dge <- DGEList(counts = count_matrix_valid, group = group_labels)

# Filter out lowly expressed genes
keep <- filterByExpr(dge, group = group_labels)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization
dge <- calcNormFactors(dge)

#---------------------------------------------------------------
# Design Matrix and voom Transformation
#---------------------------------------------------------------

design <- model.matrix(~0 + group_labels)
colnames(design) <- levels(group_labels)  # "Bad" and "Good"

v <- voom(dge, design, plot = TRUE)  # plot=TRUE shows mean-variance trend

#---------------------------------------------------------------
# Fit Linear Model and Make Contrast (Bad vs Good)
#---------------------------------------------------------------

fit <- lmFit(v, design)

# Define the comparison of interest: Bad - Good
contrast_matrix <- makeContrasts(Bad_vs_Good = Bad - Good, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

#---------------------------------------------------------------
# Extract Differentially Expressed Genes
#---------------------------------------------------------------

deg_results <- topTable(fit2, coef = "Bad_vs_Good", number = Inf, adjust.method = "BH")

# Add gene identifiers
deg_results$Gene <- rownames(deg_results)
deg_results <- deg_results[, c("Gene", setdiff(colnames(deg_results), "Gene"))]

#---------------------------------------------------------------
# Filter Significant DEGs
#---------------------------------------------------------------

deg_sig <- subset(deg_results, abs(logFC) > 0.2 & P.Value < 0.01) 
nrow(deg_sig)

#---------------------------------------------------------------
# Prepare clinical annotation for heatmap (corrected, preserve order)
#---------------------------------------------------------------

# Ensure consistent ordering
tnbc_df <- tnbc_df[order(tnbc_df$bcr_patient_barcode), ]

annotation_df <- tnbc_df %>%
  select(
    bcr_patient_barcode,
    age_at_diagnosis,
    ajcc_tumor_pathologic_pt,
    ajcc_nodes_pathologic_pn,
    ajcc_pathologic_tumor_stage,
    bad_prognosis
  )

# Rename columns
colnames(annotation_df) <- c("Patient", "Age", "Tumor", "Nodes", "Stage", "Prognosis")

# Convert prognosis logicals to categorical labels
annotation_df$Prognosis <- ifelse(annotation_df$Prognosis, "Bad", "Good")

# Convert variables to factors
annotation_df$Tumor <- as.factor(annotation_df$Tumor)
annotation_df$Nodes <- as.factor(annotation_df$Nodes)
annotation_df$Stage <- as.factor(annotation_df$Stage)
annotation_df$Prognosis <- as.factor(annotation_df$Prognosis)

# Set rownames as patient IDs (preserve order from tnbc_df)
rownames(annotation_df) <- annotation_df$Patient

# Show the cleaned dataframe
annotation_df

#---------------------------------------------------------------
# Match order of clinical data to expression matrix columns
#---------------------------------------------------------------

# Extract patient IDs from expression matrix
expr_samples <- substr(colnames(count_matrix_valid), 1, 12)

# Subset and reorder annotation to match expression columns
annotation_df <- annotation_df[match(expr_samples, annotation_df$Patient), ]

# Set rownames of annotation to patient/sample IDs
rownames(annotation_df) <- annotation_df$Patient

#---------------------------------------------------------------
# Subset expression matrix for significant DEGs
#---------------------------------------------------------------
genes_for_heatmap <- intersect(rownames(v$E), rownames(deg_sig[order(deg_sig$P.Value), ]))
mat_all <- v$E[genes_for_heatmap, ]

#---------------------------------------------------------------
# Prepare clinical annotation for heatmap with categories
#---------------------------------------------------------------

# Copy original annotation and preserve Patient IDs
annotation_clean <- annotation_df

#---------------------------------------------------------------
# Age: numeric bins (10-year increments)
#---------------------------------------------------------------

annotation_clean$Age_num <- suppressWarnings(as.numeric(as.character(annotation_clean$Age)))

annotation_clean$Age_cat <- cut(
  annotation_clean$Age_num,
  breaks = seq(20, 90, by = 10),
  include.lowest = TRUE,
  right = FALSE,
  labels = paste(seq(20, 80, by = 10), seq(29, 89, by = 10), sep = "-")
)

#---------------------------------------------------------------
# Tumor: T1-T2 vs T3-T4
#---------------------------------------------------------------

annotation_clean$Tumor_cat <- ifelse(
  grepl("^T[12]", annotation_clean$Tumor, ignore.case = TRUE), "T1-T2",
  ifelse(grepl("^T[34]", annotation_clean$Tumor, ignore.case = TRUE), "T3-T4", "Unknown")
)

#---------------------------------------------------------------
# Node: N0-N1 vs N2-N3
#---------------------------------------------------------------

annotation_clean$Node_cat <- ifelse(
  grepl("^N[01]", annotation_clean$Nodes, ignore.case = TRUE), "N0-N1",
  ifelse(grepl("^N[23]", annotation_clean$Nodes, ignore.case = TRUE), "N2-N3", "Unknown")
)

#---------------------------------------------------------------
# Stage: Stage I-II vs Stage III-IV
#---------------------------------------------------------------

annotation_clean$Stage_cat <- ifelse(
  grepl("Stage\\s*[I]{1,2}", annotation_clean$Stage, ignore.case = TRUE), "Stage I-II",
  ifelse(grepl("Stage\\s*[I]{3,4}", annotation_clean$Stage, ignore.case = TRUE), "Stage III-IV", "Unknown")
)

#---------------------------------------------------------------
# Prognosis: already categorical
#---------------------------------------------------------------

annotation_clean$Prognosis_cat <- ifelse(annotation_clean$Prognosis == "Good", "Good", "Bad")

#---------------------------------------------------------------
# Keep only categorical columns for heatmap
#---------------------------------------------------------------

heatmap_annotation <- annotation_clean[, c("Age_cat", "Tumor_cat", "Node_cat", "Stage_cat", "Prognosis_cat")]

# Convert all to factors
heatmap_annotation[] <- lapply(heatmap_annotation, as.factor)

# Set rownames to patient IDs (matching expression matrix order)
expr_samples <- substr(colnames(count_matrix_valid), 1, 12)
heatmap_annotation <- heatmap_annotation[match(expr_samples, annotation_clean$Patient), ]
rownames(heatmap_annotation) <- expr_samples

#---------------------------------------------------------------
# Prepare annotation for heatmap: Age, Tumor, Node, Stage, Prognosis
#---------------------------------------------------------------

annotation_for_heatmap <- annotation_df

# Age: bin in 10-year increments
annotation_for_heatmap$Age_group <- cut(
  as.numeric(annotation_for_heatmap$Age),
  breaks = seq(20, 90, by = 10),
  right = FALSE,
  include.lowest = TRUE
)
annotation_for_heatmap$Age_group <- factor(annotation_for_heatmap$Age_group)

# Tumor: T1-T2 vs T3-T4
annotation_for_heatmap$Tumor_group <- ifelse(
  grepl("^T[12]", annotation_for_heatmap$Tumor, ignore.case = TRUE), "T1-T2", "T3-T4"
)
annotation_for_heatmap$Tumor_group <- factor(annotation_for_heatmap$Tumor_group)

# Node: N0-N1 vs N2-N3
annotation_for_heatmap$Node_group <- ifelse(
  grepl("^N[01]", annotation_for_heatmap$Nodes, ignore.case = TRUE), "N0-N1", "N2-N3"
)
annotation_for_heatmap$Node_group <- factor(annotation_for_heatmap$Node_group)

# Stage: Stage I-II vs Stage III-IV
# Classify Stage groups including Unknown
annotation_for_heatmap$Stage_group <- ifelse(
  grepl("^Stage\\s*I{1,2}[AB]?$", annotation_for_heatmap$Stage, ignore.case = TRUE), 
  "I-II", 
  ifelse(
    grepl("^Stage\\s*III|IV", annotation_for_heatmap$Stage, ignore.case = TRUE), 
    "III-IV", 
    "Unknown"  # assign "Unknown" to anything else or missing
  )
)

annotation_for_heatmap$Stage_group <- factor(annotation_for_heatmap$Stage_group, levels = c("I-II", "III-IV", "Unknown"))

table(annotation_for_heatmap$Stage_group)

# Prognosis: already Good/Bad
annotation_for_heatmap$Prognosis <- factor(annotation_for_heatmap$Prognosis, levels = c("Good", "Bad"))

# Subset only the grouped annotation columns for heatmap
heatmap_annotation <- annotation_for_heatmap[, c("Age_group", "Tumor_group", "Node_group", "Stage_group", "Prognosis")]

# Rename columns for clarity if you want (optional)
colnames(heatmap_annotation) <- c("Age", "Tumor", "Node", "Stage", "Prognosis")

# Ensure rownames match expression matrix columns
rownames(heatmap_annotation) <- colnames(mat_all)

#---------------------------------------------------------------
# Define colors for each annotation (named)
#---------------------------------------------------------------

# Age colors: named by factor levels
age_levels <- levels(heatmap_annotation$Age)
age_colors <- colorRampPalette(c("#ffffcc", "#800026"))(length(age_levels))
names(age_colors) <- age_levels


annotation_colors <- list(
  Age       = age_colors,  # keep your existing age palette
  Tumor     = c("T1-T2" = "#A8E6CF",  
                "T3-T4" = "#1E90FF"),
  Node      = c("N0-N1" = "#D0B0FF",
                "N2-N3" = "#FFB6C1"), 
  Stage     = c("I-II" = "#C7F464",   
                "III-IV" = "#FFFACD",  
                "Unknown" = "#A6A6A6"),
  Prognosis = c("Good" = "#FA8072",    
                "Bad" = "#00CED1")     
)

#---------------------------------------------------------------
# Plot heatmap
#---------------------------------------------------------------
pheatmap(
  mat_all,
  annotation_col = heatmap_annotation,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  width = 12,   # adjust width (in inches)
  height = 8,    # adjust height
  main = "DEGs: Bad vs Good Prognosis"
)

#---------------------------------------------------------------
# Extract Data for Integration
#---------------------------------------------------------------

# Use your list of matched patients (already curated)
match_patients <- final_patient_ids

# Extract patient columns matching your curated list
sample_patients <- substr(colnames(v$E), 1, 12)  # first 12 chars of sample barcodes
selected_cols <- which(sample_patients %in% match_patients)

# Extract rows for DEGs (use Ensembl IDs with versions)
deg_ensembl_ids <- deg_sig$Gene  # still contains ENSG IDs with versions
selected_rows <- which(rownames(v$E) %in% deg_ensembl_ids)

# Subset the expression matrix
expr_deg <- v$E[selected_rows, selected_cols]

# Get gene annotation from the SummarizedExperiment object
gene_info <- rowData(tcga_brca_transcriptomic)

# Add gene symbols for readability
gene_symbols <- gene_info$gene_name[match(rownames(expr_deg), gene_info$gene_id)]
rownames(expr_deg) <- gene_symbols

# Check dimensions
dim(expr_deg)  # rows = DEGs, columns = selected patients
head(expr_deg)

# Make sure rownames (genes) and colnames (patients) are set
# expr_deg: rows = genes, columns = patient IDs
expr_deg_df <- as.data.frame(expr_deg)
expr_deg_df$Gene <- rownames(expr_deg_df)  # move gene names to a column

# Pivot longer to get long format (Patient × Gene × log2Expr)
expr_long <- expr_deg_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Patient_ID",
    values_to = "Log2CPM"  # voom-normalized expression
  )

# Standardize TCGA identifiers to 12-character patient-level IDs
expr_long$Patient_ID <- substr(expr_long$Patient_ID, 1, 12)

# Save to Excel
write_xlsx(expr_long, "TNBC_For.Integration_Transcriptomics.xlsx")

