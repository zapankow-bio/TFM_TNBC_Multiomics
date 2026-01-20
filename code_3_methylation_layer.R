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

# Output to excel (optional)
#tnbc_patient_ids_df <- as.data.frame(tnbc_patient_ids)
#write_xlsx(tnbc_patient_ids_df, "TNBC_Clinical_Patients.xlsx")

# Print how many triple negative patients were found
cat("Number of triple negative breast cancer patients:", length(tnbc_patient_ids), "\n")

# View patient IDs
print(tnbc_patient_ids)

# Result: 116 patients have TNBC.

#-----------------------------------------------
# DNA Methylation Layer
#-----------------------------------------------

# Create a query to access DNA Methylation data from TCGA-BRCA project
query_methylation <- GDCquery(
  project= "TCGA-BRCA", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  sample.type = "Primary Tumor"
)

# Download the DNA Methylation data files based on the query above
GDCdownload(query_methylation, method = "api", files.per.chunk = 10)

# Prepare and load the downloaded DNA Methylation data into R as a SummarizedExperiment object
tcga_brca_methylation <- GDCprepare(query_methylation, summarizedExperiment = TRUE)

# Extract methylation beta value matrix (CpGs x samples)
beta_values <- assay(tcga_brca_methylation)

# Extract full sample barcodes from the methylation data
methylation_sample_barcodes <- colnames(beta_values)

# Extract patient barcodes (first 12 characters)
methylation_patient_barcodes <- substr(methylation_sample_barcodes, 1, 12)

# Subset: Only those patient barcodes that are in the TNBC list
tnbc_methylation_match <- methylation_patient_barcodes %in% tnbc_patient_ids

# Get the sample barcodes and corresponding patient barcodes
tnbc_methylation_sample_barcodes <- methylation_sample_barcodes[tnbc_methylation_match]
tnbc_methylation_patient_barcodes <- methylation_patient_barcodes[tnbc_methylation_match]

# Create a data frame of samples + patients
tnbc_methylation_df <- data.frame(
  Sample_Barcode = tnbc_methylation_sample_barcodes,
  Patient_Barcode = tnbc_methylation_patient_barcodes,
  stringsAsFactors = FALSE
)

# Find duplicate patients
dup_patients <- tnbc_methylation_df$Patient_Barcode[duplicated(tnbc_methylation_df$Patient_Barcode)]
unique(dup_patients)

#-----------------------------------------------
# Resolve Duplicates by Comparing Quality
#-----------------------------------------------

# This assumes beta_values contains all methylation samples
colnames(beta_values) <- methylation_sample_barcodes

# Store selected best samples
best_samples <- c()

# For each duplicated patient
for (pid in unique(dup_patients)) {
  
  # Get all sample barcodes for this patient
  sample_barcodes <- tnbc_methylation_df$Sample_Barcode[tnbc_methylation_df$Patient_Barcode == pid]
  
  # Subset beta values
  beta_subset <- beta_values[, sample_barcodes, drop = FALSE]
  
  # Calculate non-missing CpG counts per sample
  non_na_counts <- colSums(!is.na(beta_subset))
  
  # Pick the sample with the most non-NA CpGs
  best_sample <- names(which.max(non_na_counts))
  
  best_samples <- c(best_samples, best_sample)
}

# Get patients who had only one sample
unique_patients <- tnbc_methylation_df$Patient_Barcode[!duplicated(tnbc_methylation_df$Patient_Barcode) & 
                                                         !tnbc_methylation_df$Patient_Barcode %in% dup_patients]

# Get their corresponding sample barcodes
unique_samples <- tnbc_methylation_df$Sample_Barcode[tnbc_methylation_df$Patient_Barcode %in% unique_patients]

# Combine all selected samples
final_selected_samples <- c(best_samples, unique_samples)

# Subset beta matrix to selected samples
beta_values_tnbc_filtered <- beta_values[, final_selected_samples]

# Extract full sample barcodes from the methylation data
methylation_sample_barcodes_unique <- colnames(beta_values_tnbc_filtered)

# Extract patient barcodes (first 12 characters)
methylation_patient_barcodes_unique <- substr(methylation_sample_barcodes_unique, 1, 12)

# Subset: Only those patient barcodes that are in the TNBC list
tnbc_methylation_match_unique <- methylation_patient_barcodes_unique %in% tnbc_patient_ids

# Get the sample barcodes and corresponding patient barcodes
tnbc_methylation_sample_barcodes_final <- methylation_sample_barcodes_unique[tnbc_methylation_match_unique]
tnbc_methylation_patient_barcodes_final <- methylation_patient_barcodes_unique[tnbc_methylation_match_unique]

# Result: 83 TNBC patients have DNA Methylation data

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
# Binary DNA Methylation Matrix and Heatmap
#---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

#---------------------------------------------------------------
# 1. Match methylation samples to curated 60 TNBC patients
#---------------------------------------------------------------

methyl_df <- data.frame(
  Sample = colnames(beta_values_tnbc_filtered),
  Patient_ID = substr(colnames(beta_values_tnbc_filtered), 1, 12),
  stringsAsFactors = FALSE
)

# Keep only curated 60 patients
keep_idx <- methyl_df$Patient_ID %in% final_patient_ids
beta_60 <- beta_values_tnbc_filtered[, keep_idx, drop = FALSE]
methyl_df <- methyl_df[keep_idx, ]

# Order columns by patient barcode
beta_60 <- beta_60[, order(methyl_df$Patient_ID)]
methyl_df <- methyl_df[order(methyl_df$Patient_ID), ]

#---------------------------------------------------------------
# 2. Create binary methylation matrix
#---------------------------------------------------------------
# Rule: 1 = hyper OR hypo; 0 = normal range
binary_mat <- ifelse(beta_60 < 0.2 | beta_60 > 0.8, 1, 0)

#---------------------------------------------------------------
# 3. Filter CpGs to reduce noise
#---------------------------------------------------------------

min_samples <- 5
keep_cpgs <- rowSums(binary_mat) >= min_samples
binary_mat_filt <- binary_mat[keep_cpgs, , drop = FALSE]

cat("CpGs retained:", nrow(binary_mat_filt), "\n")

#---------------------------------------------------------------
# 4. Add prognosis annotation
#---------------------------------------------------------------

annot <- data.frame(
  Patient_ID = methyl_df$Patient_ID,
  Sample = methyl_df$Sample,
  stringsAsFactors = FALSE
)

annot <- annot %>%
  left_join(
    tnbc_df %>% dplyr::select(bcr_patient_barcode, bad_prognosis) %>%
      rename(Patient_ID = bcr_patient_barcode),
    by = "Patient_ID"
  ) %>%
  mutate(Prognosis = ifelse(bad_prognosis, "Bad", "Good"))

annotation_col <- data.frame(Prognosis = annot$Prognosis)
rownames(annotation_col) <- annot$Sample

#---------------------------------------------------------------
# 5. Map CpGs to genes
#---------------------------------------------------------------
# Make sure you have the annotation loaded previously
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Only keep CpGs present in filtered binary matrix
anno_sub <- anno[rownames(binary_mat_filt), ]
genes <- sapply(strsplit(as.character(anno_sub$UCSC_RefGene_Name), ";"), `[`, 1)

binary_mat_sig_num <- binary_mat_filt

#---------------------------------------------------------------
# 6. Aggregate CpGs to gene-level (1 if any CpG abnormal)
#---------------------------------------------------------------

split_list <- split(as.data.frame(binary_mat_sig_num), genes)

gene_mat <- do.call(rbind, lapply(split_list, function(df_gene) {
  apply(df_gene, 2, function(x) as.numeric(max(x, na.rm = TRUE)))
}))

rownames(gene_mat) <- names(split_list)

#---------------------------------------------------------------
# Pipeline for Differentially Methylation Analysis
#---------------------------------------------------------------

#---------------------------------------------------------------
# 1. Load libraries
#---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

#---------------------------------------------------------------
# 2. Prepare Beta Matrix for Curated 60 TNBC Patients
#---------------------------------------------------------------

methyl_df <- data.frame(
  Sample = colnames(beta_values_tnbc_filtered),
  Patient_ID = substr(colnames(beta_values_tnbc_filtered), 1, 12),
  stringsAsFactors = FALSE
)

keep_idx <- methyl_df$Patient_ID %in% final_patient_ids
beta_60 <- beta_values_tnbc_filtered[, keep_idx, drop = FALSE]
methyl_df <- methyl_df[keep_idx, ]

# Order columns by patient barcode
beta_60 <- beta_60[, order(methyl_df$Patient_ID)]
methyl_df <- methyl_df[order(methyl_df$Patient_ID), ]

#---------------------------------------------------------------
# 3. Match Clinical Prognosis
#---------------------------------------------------------------

tnbc_matched_df <- tnbc_df[match(methyl_df$Patient_ID, tnbc_df$bcr_patient_barcode), ]
stopifnot(all(methyl_df$Patient_ID == tnbc_matched_df$bcr_patient_barcode))

group <- factor(ifelse(tnbc_matched_df$bad_prognosis, "Bad", "Good"))

#---------------------------------------------------------------
# 4. Convert Beta → M-values
#---------------------------------------------------------------

M_values <- log2(beta_60 / (1 - beta_60))
M_values[!is.finite(M_values)] <- NA

# Keep CpGs present in ≥80% of samples
keep_cpgs <- rowMeans(!is.na(M_values)) >= 0.8
M_values_filtered <- M_values[keep_cpgs, ]

#---------------------------------------------------------------
# 5. Differential Methylation (Bad vs Good)
#---------------------------------------------------------------

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(M_values_filtered, design)
contrast.matrix <- makeContrasts(Bad_vs_Good = Bad - Good, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Δβ for interpretability
mean_beta <- t(apply(beta_60[rownames(fit2), ], 1, function(x) tapply(x, group, mean, na.rm = TRUE)))
topDMRs <- topTable(fit2, coef = "Bad_vs_Good", number = Inf, adjust.method = "BH")
topDMRs$deltaBeta <- mean_beta[, "Bad"] - mean_beta[, "Good"]

#---------------------------------------------------------------
# 6. Annotate CpGs to Genes
#---------------------------------------------------------------

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotated <- merge(topDMRs, anno, by.x = "row.names", by.y = "Name", all.x = TRUE)

# Clean Gene column
annotated$Gene <- sapply(strsplit(annotated$UCSC_RefGene_Name, ";", fixed = TRUE), `[`, 1)
annotated$Gene[annotated$Gene == "" | is.na(annotated$Gene)] <- NA

# Convert DFrame/list columns → data.frame with numeric vectors
annotated_clean <- annotated %>%
  as.data.frame() %>%
  filter(!is.na(Gene))

annotated_clean$logFC <- as.numeric(annotated_clean$logFC)
annotated_clean$deltaBeta <- as.numeric(annotated_clean$deltaBeta)
annotated_clean$P.Value <- as.numeric(annotated_clean$P.Value)

# Restrict CpGs before collapsing to genes
promoter_regions <- c("TSS200", "TSS1500", "5'UTR", "1stExon")

annotated_promoter <- annotated_clean %>%
  filter(UCSC_RefGene_Group %in% promoter_regions)

# Run gene-level aggregation using only promoter CpGs
diff_methylation_promoter <- annotated_promoter %>%
  group_by(Gene) %>%
  summarise(
    logFC     = mean(logFC, na.rm = TRUE),
    deltaBeta = mean(deltaBeta, na.rm = TRUE),
    P.Value   = min(P.Value, na.rm = TRUE),
    nCpGs     = n(),
    .groups = "drop"
  )

#---------------------------------------------------------------
# Top Hypermethylated and Hypomethylated Genes
#---------------------------------------------------------------

hyper_dmgs_promoter <- diff_methylation_promoter %>%
  filter(logFC > 0.5, P.Value < 0.05)

hypo_dmgs_promoter <- diff_methylation_promoter %>%
  filter(logFC < -0.5, P.Value < 0.05)

n_hyper <- nrow(hyper_dmgs_promoter)
n_hypo  <- nrow(hypo_dmgs_promoter)

cat("Hypermethylated genes:", n_hyper, "\n")
cat("Hypomethylated genes:", n_hypo, "\n")
cat("Total differentially methylated genes:", n_hyper + n_hypo, "\n")

#---------------------------------------------------------------
# Heatmap for Top Hypermethylated and Hypomethylated Genes
#---------------------------------------------------------------

library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

#---------------------------------------------------------------
# 1. Select all differentially methylated genes
#---------------------------------------------------------------

dm_genes <- diff_methylation_promoter %>%
  filter(P.Value < 0.05 & abs(logFC) > 0.5) %>%
  pull(Gene) %>%
  unique()

length(dm_genes)  # number of DMGs in heatmap

#---------------------------------------------------------------
# 2. Select all differentially methylated genes
#---------------------------------------------------------------

cpg_gene_map <- data.frame(
  CpG  = rownames(topDMRs),
  Gene = sapply(strsplit(anno[rownames(topDMRs), "UCSC_RefGene_Name"], ";"), `[`, 1),
  stringsAsFactors = FALSE
)

# Keep only DM genes
cpg_gene_map <- cpg_gene_map %>%
  filter(Gene %in% dm_genes)

# Keep CpGs present in beta matrix
cpg_gene_map <- cpg_gene_map %>%
  filter(CpG %in% rownames(beta_60))

#---------------------------------------------------------------
# 3. Collapse CpGs → gene-level mean β
#---------------------------------------------------------------

beta_long <- beta_60[cpg_gene_map$CpG, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  left_join(cpg_gene_map, by = "CpG")

beta_gene <- beta_long %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

heatmap_mat <- beta_gene %>%
  column_to_rownames("Gene") %>%
  as.matrix()

#---------------------------------------------------------------
# 4. Prognosis annotation (column-aligned)
#---------------------------------------------------------------

group_vec <- tnbc_matched_df$bad_prognosis
names(group_vec) <- methyl_df$Sample
group_vec <- group_vec[colnames(heatmap_mat)]

group_factor <- factor(ifelse(group_vec, "Bad", "Good"),
                       levels = c("Good", "Bad"))

ha <- HeatmapAnnotation(
  Prognosis = group_factor,
  col = list(Prognosis = c("Good" = "#FA8072", "Bad" = "#00CED1"))
)

#---------------------------------------------------------------
# 5. Draw heatmap (all DMGs)
#---------------------------------------------------------------

ht_all_dmgs <- Heatmap(
  heatmap_mat,
  name = "Beta",
  col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  row_names_gp = gpar(fontsize = 6)
)

draw(ht_all_dmgs)

#----------------------------------------------------------------
# Extract Data for Integration (Top Hyper and Hypomethylated genes)
#----------------------------------------------------------------

# Top hyper- and hypomethylated genes
top_genes <- unique(c(hyper_dmgs_promoter$Gene, hypo_dmgs_promoter$Gene))

# Map CpGs to these top genes
cpg_top_map <- annotated_clean %>%
  select(CpG = Row.names, Gene) %>%
  filter(Gene %in% top_genes) %>%
  filter(CpG %in% rownames(beta_60))

# Match curated patients
match_patients <- final_patient_ids
sample_patients <- substr(colnames(beta_60), 1, 12)
selected_cols <- which(sample_patients %in% match_patients)

# Subset beta matrix (top CpGs × selected patients)
beta_sub <- beta_60[rownames(beta_60) %in% cpg_top_map$CpG, selected_cols, drop = FALSE]

# Collapse CpGs → gene-level mean beta
beta_df <- as.data.frame(beta_sub)
beta_df$CpG <- rownames(beta_df)
beta_df <- beta_df %>%
  left_join(cpg_top_map, by = "CpG") %>%
  select(-CpG) %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Wide matrix
beta_wide <- beta_df %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Long format
beta_long <- beta_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Patient_ID",
    values_to = "Beta"
  )

# Standardize patient IDs to 12-character format
beta_long$Patient_ID <- factor(substr(beta_long$Patient_ID, 1, 12), levels = match_patients)
beta_long <- beta_long %>% arrange(Patient_ID, Gene)

# Optional export
write_xlsx(beta_long, "TNBC_For.Integration_Methylation.xlsx")
