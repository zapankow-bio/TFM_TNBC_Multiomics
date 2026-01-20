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
GDCdownload(query_clinical)

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
# After analysis of all the layers and clinical data, the number of TNBC patients
# ...that had data for all the layers was 60.

#---------------------------------------------------------------
# Load the 60 TNBC patient list from Excel
#---------------------------------------------------------------
tnbc_final <- read_xlsx("TNBC_Patients_Final.xlsx", col_names = FALSE)

# Ensure patient barcodes are uppercase for consistency
final_patient_ids <- toupper(tnbc_final[[1]])  # Assuming first column contains patient IDs

#---------------------------------------------------------------
# Curated TNBC Patients
#---------------------------------------------------------------

# Select only the 60 curated TNBC patients
tnbc <- patient_data[patient_data$bcr_patient_barcode %in% final_patient_ids, ]

# Convert to a data.frame
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
# Extract Data for Integration
#---------------------------------------------------------------

# Create patient-level prognosis table
patient_prognosis <- tnbc_df %>%
  dplyr::select(bcr_patient_barcode, bad_prognosis) %>%
  dplyr::mutate(
    Patient_ID = bcr_patient_barcode,
    Prognosis = ifelse(bad_prognosis, "Bad", "Good")
  ) %>%
  dplyr::select(Patient_ID, Prognosis)

# Save to Excel
write_xlsx(patient_prognosis, "TNBC_For.Integration_Clinical.xlsx")


