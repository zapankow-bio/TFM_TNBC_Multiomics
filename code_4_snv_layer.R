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
# Simple Nucleotide Variation Layer
#-----------------------------------------------

# Create a query to access SNV data from TCGA-BRCA project
query_SNV <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

# Download the SNV data files based on the query above
GDCdownload(query_SNV, method = "api", files.per.chunk = 10)

# Prepare and load the downloaded SNV data into R as a SummarizedExperiment object
tcga_brca_SNV <- GDCprepare(query_SNV, summarizedExperiment = TRUE)

# Extract SNV Sample Barcodes
SNV_barcodes <- tcga_brca_SNV$Tumor_Sample_Barcode

# Extract first 12 characters (patient IDs)
SNV_patient_ids <- substr(SNV_barcodes, 1, 12)

# Add patient ID column
tcga_brca_SNV$Patient_ID <- SNV_patient_ids

# Filter SNVs for TNBC patients
tcga_brca_SNV_tnbc <- subset(tcga_brca_SNV, Patient_ID %in% tnbc_patient_ids)

# Check how many mutation entries you have for TNBCs
cat("Number of SNV entries for TNBC patients:", nrow(tcga_brca_SNV_tnbc), "\n")

# Optional: check number of unique TNBC patients with SNV data
length(unique(tcga_brca_SNV_tnbc$Patient_ID))

#-----------------------------------------------
# Resolve Duplicates by Comparing Quality
#-----------------------------------------------

# Create a data frame from the SNV data
snv_df <- as.data.frame(tcga_brca_SNV_tnbc)

# Remove entries with missing gene symbols
snv_df <- snv_df[!is.na(snv_df$Hugo_Symbol), ]

# Filter for non-silent mutations
non_silent <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", 
                "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site")
snv_df <- snv_df[snv_df$Variant_Classification %in% non_silent, ]

# Count mutations per Tumor_Sample_Barcode
mutation_counts <- as.data.frame(table(snv_df$Tumor_Sample_Barcode))
colnames(mutation_counts) <- c("Tumor_Sample_Barcode", "Mutation_Count")

# Add Patient_IDs
mutation_counts$Patient_ID <- substr(mutation_counts$Tumor_Sample_Barcode, 1, 12)

# For each patient, keep the sample with the highest mutation count
library(dplyr)
best_samples <- mutation_counts %>%
  group_by(Patient_ID) %>%
  arrange(desc(Mutation_Count)) %>%
  slice(1) %>%
  ungroup()

# Filter original SNV data to keep only best samples
tcga_brca_SNV_tnbc_filtered <- snv_df[snv_df$Tumor_Sample_Barcode %in% best_samples$Tumor_Sample_Barcode, ]

# Optional: confirm dimensions
cat("Number of mutation entries after filtering:", nrow(tcga_brca_SNV_tnbc_filtered), "\n")
cat("Number of unique patients:", length(unique(tcga_brca_SNV_tnbc_filtered$Patient_ID)), "\n")
cat("Number of unique samples:", length(unique(tcga_brca_SNV_tnbc_filtered$Tumor_Sample_Barcode)), "\n")

# Result: 98 TNBC patients have SNV data.

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

#---------------------------------------------------------------
# Load curated 60 patients
#---------------------------------------------------------------

# Count mutations per sample
mut_burden <- tcga_brca_SNV_tnbc_filtered %>%
  group_by(Tumor_Sample_Barcode, Patient_ID) %>%
  summarise(Mutation_Count = n(), .groups = "drop")

# Add group (Good vs Bad prognosis) by merging with tnbc_df
mut_burden <- mut_burden %>%
  left_join(
    tnbc_df %>%
      select(bcr_patient_barcode, bad_prognosis) %>%
      rename(Patient_ID = bcr_patient_barcode),
    by = "Patient_ID"
  ) %>%
  mutate(
    Group = ifelse(bad_prognosis, "Bad", "Good")
  )

# final_patient_ids contains your curated 60 TNBC patient IDs
mut_burden_60 <- mut_burden %>%
  filter(substr(Tumor_Sample_Barcode, 1, 12) %in% final_patient_ids)

# Check how many patients per prognosis group in curated 60
mut_burden_60 %>%
  group_by(Group) %>%
  summarise(num_patients = n())

#---------------------------------------------------------------
# Binary SNV Heatmap with Prognosis Annotation
#---------------------------------------------------------------

# Libraries
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Prepare mutation matrix for curated 60 samples

# Ensure patient IDs are first 12 chars
tcga_brca_SNV_tnbc_filtered$Patient_ID <- substr(tcga_brca_SNV_tnbc_filtered$Tumor_Sample_Barcode, 1, 12)

# Keep only curated 60 patients
filtered_samples <- tcga_brca_SNV_tnbc_filtered[
  substr(tcga_brca_SNV_tnbc_filtered$Tumor_Sample_Barcode, 1, 12) %in% final_patient_ids, ]

# Create binary mutation presence table: rows = Hugo_Symbol, cols = Tumor_Sample_Barcode
mut_mat <- filtered_samples %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Tumor_Sample_Barcode,
              values_from = present,
              values_fill = list(present = 0))

# Convert to matrix
genes <- mut_mat$Hugo_Symbol
mut_mat_matrix <- as.matrix(mut_mat[, -1, drop = FALSE])
rownames(mut_mat_matrix) <- genes

# Filter genes mutated in at least min_samples patients
min_samples <- 5
keep_genes <- rowSums(mut_mat_matrix) >= min_samples
mut_mat_matrix_filt <- mut_mat_matrix[keep_genes, , drop = FALSE]

# Create column annotation (Prognosis)
sample_patient_map <- data.frame(
  Tumor_Sample_Barcode = colnames(mut_mat_matrix_filt),
  Patient_ID = substr(colnames(mut_mat_matrix_filt), 1, 12),
  stringsAsFactors = FALSE
)

sample_patient_map <- sample_patient_map %>%
  left_join(
    tnbc_df %>% dplyr::select(bcr_patient_barcode, bad_prognosis) %>%
      dplyr::rename(Patient_ID = bcr_patient_barcode),
    by = "Patient_ID"
  ) %>%
  mutate(Prognosis = ifelse(is.na(bad_prognosis), "Unknown",
                            ifelse(bad_prognosis, "Bad", "Good")))

annotation_col <- data.frame(Prognosis = sample_patient_map$Prognosis)
rownames(annotation_col) <- sample_patient_map$Tumor_Sample_Barcode

# Ensure column order matches annotation
mut_mat_matrix_filt <- mut_mat_matrix_filt[, rownames(annotation_col), drop = FALSE]

# Binary matrix for heatmap
mat_bin <- ifelse(mut_mat_matrix_filt > 0, 1, 0)

# Define binary colors
bin_colors <- c("0" = "white", "1" = "red")


# Heatmap
ht <- Heatmap(
  mat_bin,
  name = "Non-Syn. Mutation",
  col = bin_colors,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  
  # Row names bold
  row_names_gp = gpar(fontsize = 6, fontface = "bold"),
  
  # Cell borders (grid lines)
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(
      x, y, width, height,
      gp = gpar(col = "gray70", fill = fill, lwd = 0.5)
    )
  },
  
  # Top annotation for Prognosis
  top_annotation = HeatmapAnnotation(
    df = annotation_col,
    col = list(
      Prognosis = c("Good" = "#FA8072", "Bad" = "#00CED1", "Unknown" = "gray80")
    )
  ),
  
  # Legend customization
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("No mutation", "Mutation"),
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

# Draw heatmap
draw(ht)

#---------------------------------------------------------------
# Extract Data for Integration
#---------------------------------------------------------------

# Start with filtered SNV data (best sample per patient)
snv_int <- tcga_brca_SNV_tnbc_filtered %>%
  mutate(Patient_ID = substr(Tumor_Sample_Barcode, 1, 12)) %>%
  filter(Patient_ID %in% final_patient_ids) %>%
  dplyr::select(Patient_ID, Hugo_Symbol, Variant_Classification) %>%
  distinct()

# Get all unique genes and patients
all_genes <- unique(snv_int$Hugo_Symbol)
all_patients <- unique(snv_int$Patient_ID)

# Create full patient × gene grid
full_grid <- expand.grid(
  Patient_ID = all_patients,
  Gene = all_genes,
  stringsAsFactors = FALSE
)

# Merge with mutation data
snv_long_full <- full_grid %>%
  left_join(
    snv_int %>% select(Patient_ID, Gene = Hugo_Symbol) %>% mutate(Mutation = 1),
    by = c("Patient_ID", "Gene")
  ) %>%
  mutate(Mutation = ifelse(is.na(Mutation), 0, Mutation))  # non-mutated = 0

# Rearrange columns
snv_long_full <- snv_long_full %>%
  select(Gene, Patient_ID, Mutation)

# Count number of patients mutated per gene
gene_mut_counts <- snv_long_full %>%
  group_by(Gene) %>%
  summarise(Patient_Count = sum(Mutation), .groups = "drop")

# Keep genes mutated in at least 3 patients
genes_keep <- gene_mut_counts %>%
  filter(Patient_Count >= 3) %>%
  pull(Gene)

# Filter SNV long table
snv_long_filtered <- snv_long_full %>%
  filter(Gene %in% genes_keep)

# Save to Excel
write_xlsx(snv_long_filtered, "TNBC_For.Integration_SNV.xlsx")

