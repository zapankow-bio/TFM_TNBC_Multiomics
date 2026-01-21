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
library(ComplexHeatmap)
library(tidyr)
library(dplyr)

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
#GDCdownload(query_clinical)

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

# Print how many triple negative patients were found
cat("Number of triple negative breast cancer patients:", length(tnbc_patient_ids), "\n")

# View patient IDs
print(tnbc_patient_ids)

# Result: 116 patients have TNBC.

#-----------------------------------------------
# Copy Number Variation Layer
#-----------------------------------------------

# Create a query to access CNV data from TCGA-BRCA project
query_CNV <- GDCquery(
  project = "TCGA-BRCA",                             
  data.category = "Copy Number Variation",    
  data.type = "Gene Level Copy Number",
  access = "open",
  sample.type = "Primary Tumor"
)

# Download files if not already downloaded
GDCdownload(query_CNV, method = "api", files.per.chunk = 10)

# Prepare the CNV data
tcga_brca_CNV <- GDCprepare(query_CNV, summarizedExperiment = TRUE)

# Extract CNV sample barcodes
cnv_sample_barcodes <- colnames(assay(tcga_brca_CNV))

# Extract patient barcodes (first 12 characters)
cnv_patient_ids <- substr(cnv_sample_barcodes, 1, 12)

# Subset CNV data for TNBC patients
tnbc_cnv_match <- cnv_patient_ids %in% tnbc_patient_ids

# Filter columns (samples) to only TNBC patients
cnv_matrix_tnbc <- assay(tcga_brca_CNV)[, tnbc_cnv_match]

# Update sample and patient barcodes
tnbc_cnv_sample_barcodes <- colnames(cnv_matrix_tnbc)
tnbc_cnv_patient_ids <- substr(tnbc_cnv_sample_barcodes, 1, 12)

# Count how many CNV samples per TNBC patient
table_tnbc_cnv <- table(tnbc_cnv_patient_ids)

# Identify duplicate patients
dup_patients <- names(table_tnbc_cnv[table_tnbc_cnv > 1])
cat("Duplicated TNBC patients in CNV data:", length(dup_patients), "\n")

# Output number of TNBC patients with CNV data
cat("Number of TNBC patients with CNV data:", length(tnbc_cnv_patient_ids), "\n")

# Result: there are 106 TNBC patients with CNV data. There are no duplicates.

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

# Early death (<5 years = 1825 days)
tnbc_df$bad_death <- ifelse(!is.na(tnbc_df$death_days_to) & tnbc_df$death_days_to < 1825, TRUE, FALSE)

# Large/advanced tumor (T3 or T4, including subtypes like T4a)
tnbc_df$bad_tumor <- grepl("^T[34]", tnbc_df$ajcc_tumor_pathologic_pt, ignore.case = TRUE)

# Extensive nodal involvement (N2 or N3, including subtypes like N3a)
tnbc_df$bad_nodes <- grepl("^N[23]", tnbc_df$ajcc_nodes_pathologic_pn, ignore.case = TRUE)

# Advanced pathologic stage (Stage III or IV, including subtypes like IIIA)
tnbc_df$bad_stage <- grepl("Stage\\s*(III|IV)", tnbc_df$ajcc_pathologic_tumor_stage, ignore.case = TRUE)

# Combine tumor, nodes, and stage into one feature flag
tnbc_df$bad_prognosis <- tnbc_df$bad_death | tnbc_df$bad_tumor | tnbc_df$bad_nodes | tnbc_df$bad_stage

# Count how many patients have bad prognosis
num_bad <- sum(tnbc_df$bad_prognosis, na.rm = TRUE)
cat("Number of TNBC patients with bad prognosis:", num_bad, "\n")

# See patient IDs with bad prognosis
bad_patients <- tnbc_df$bcr_patient_barcode[tnbc_df$bad_prognosis]
print(bad_patients)

# How the CNV data is coded
# 0 → Homozygous deletion
# 1 → Heterozygous deletion
# 2 → Neutral/diploid (most common, ~1.17M counts)
# 3 → Low-level gain (trisomy)
# 4–6 → High-level gains / amplifications
# 7–16 → Very high-level amplifications (rare, only a few genes)
# Most of the CNV data is normal (2 copies) or low-level gains (3–4 copies).

#---------------------------------------------------------------
# Binary CNV Heatmap for Top 50 Genes with Prognosis Annotation
#---------------------------------------------------------------

# Libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(grid)

# Load curated CNV matrix (already prepared)
# Assume cnv_matrix_tnbc is CNV data (genes x samples) for TNBC patients
# Assume final_patient_ids contains your 60 curated TNBC patient IDs
# Assume tnbc_df contains curated 60 patients with bad_prognosis flag

# Filter CNV matrix for curated 60 patients
cnv_patient_ids <- substr(colnames(cnv_matrix_tnbc), 1, 12)
cnv_matrix_60 <- cnv_matrix_tnbc[, cnv_patient_ids %in% final_patient_ids]

# Remove rows (genes) with any NA
cnv_matrix_60_noNA <- cnv_matrix_60[!apply(cnv_matrix_60, 1, function(x) any(is.na(x))), ]

# Convert CNV to binary
cnv_matrix_binary <- cnv_matrix_60_noNA
cnv_matrix_binary[cnv_matrix_binary == 2] <- 0  # normal
cnv_matrix_binary[cnv_matrix_binary != 0] <- 1  # any gain/loss -> 1

# Filter genes altered in at least 5 patients
min_patients <- 5
gene_alter_counts <- rowSums(cnv_matrix_binary)
cnv_matrix_binary_filt <- cnv_matrix_binary[gene_alter_counts >= min_patients, ]

# Map Ensembl IDs to Hugo Gene Symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ens_ids <- gsub("\\..*$", "", rownames(cnv_matrix_binary_filt))  # remove version numbers

mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ens_ids,
  mart = ensembl
)

map_vector <- mapping$hgnc_symbol
names(map_vector) <- mapping$ensembl_gene_id

gene_symbols <- map_vector[ens_ids]

# Keep only genes with official HGNC symbols
keep_idx <- !is.na(gene_symbols) & gene_symbols != ""
cnv_matrix_binary_filt <- cnv_matrix_binary_filt[keep_idx, ]
gene_symbols <- gene_symbols[keep_idx]

# Select top 50 most variable genes
gene_var <- apply(cnv_matrix_binary_filt, 1, var)
top_genes <- order(gene_var, decreasing = TRUE)[1:min(50, length(gene_var))]
cnv_matrix_top50 <- cnv_matrix_binary_filt[top_genes, ]
rownames(cnv_matrix_top50) <- gene_symbols[top_genes]

# Column annotation: Prognosis
cnv_patient_ids_top50 <- substr(colnames(cnv_matrix_top50), 1, 12)
annotation_col <- data.frame(
  Prognosis = factor(
    ifelse(
      is.na(tnbc_df$bad_prognosis[match(cnv_patient_ids_top50, tnbc_df$bcr_patient_barcode)]),
      "Unknown",
      ifelse(
        tnbc_df$bad_prognosis[match(cnv_patient_ids_top50, tnbc_df$bcr_patient_barcode)],
        "Bad",
        "Good"
      )
    ),
    levels = c("Good", "Bad", "Unknown")
  )
)
rownames(annotation_col) <- colnames(cnv_matrix_top50)

# Define binary colors
bin_colors <- c("0" = "white", "1" = "red")

# Draw CNV heatmap
ht_cnv <- Heatmap(
  cnv_matrix_top50,
  name = "CNV Alteration",
  col = bin_colors,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  
  # Bold gene names
  row_names_gp = gpar(fontsize = 6, fontface = "bold"),
  
  # Cell borders (grid lines)
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(
      x, y, width, height,
      gp = gpar(col = "gray70", fill = fill, lwd = 0.5)
    )
  },
  
  # Top annotation for prognosis
  top_annotation = HeatmapAnnotation(
    df = annotation_col,
    col = list(
      Prognosis = c("Good" = "#FA8072", "Bad" = "#00CED1")
    )
  ),
  
  # Legend customization
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("Normal", "Abnormal"),
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

# Draw heatmap
draw(ht_cnv)

#---------------------------------------------------------------
# Extract Data for Integration
#---------------------------------------------------------------

# Remove version suffix from CNV matrix Ensembl IDs
ensembl_ids_clean <- gsub("\\..*$", "", rownames(cnv_matrix_binary_filt))

# Clean gene_symbols vector
gene_symbols_clean <- gene_symbols
names(gene_symbols_clean) <- gsub("\\..*$", "", names(gene_symbols_clean))

# Keep only genes that have a valid HGNC symbol
keep_idx <- ensembl_ids_clean %in% names(gene_symbols_clean)
cnv_matrix_binary_filt <- cnv_matrix_binary_filt[keep_idx, ]
ensembl_ids_clean <- ensembl_ids_clean[keep_idx]

# Map HGNC symbols
gene_symbols_mapped <- gene_symbols_clean[ensembl_ids_clean]

# Assign HGNC symbols as rownames
rownames(cnv_matrix_binary_filt) <- gene_symbols_mapped

# Compute variance and select top 300 most variable genes
gene_var <- apply(cnv_matrix_binary_filt, 1, var)
top_n <- 300
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(top_n, length(gene_var))]

# Subset CNV matrix for top genes
cnv_top <- cnv_matrix_binary_filt[top_genes, ]

# Subset for curated TNBC patients
match_patients <- final_patient_ids  # 12-character TCGA IDs
sample_patients <- substr(colnames(cnv_top), 1, 12)
selected_cols <- which(sample_patients %in% match_patients)

cnv_sel <- cnv_top[, selected_cols]

# Convert to long format (Gene × Patient_ID × CNV)
cnv_df <- as.data.frame(cnv_sel)
cnv_df$Gene <- rownames(cnv_df)

cnv_long <- cnv_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Patient_ID",
    values_to = "CNV"
  )

# Standardize patient IDs
cnv_long$Patient_ID <- substr(cnv_long$Patient_ID, 1, 12)

# Save to Excel
write_xlsx(cnv_long, "TNBC_For.Integration_CNV.xlsx")

