# TFM – Integración Multiómica en Cáncer de Mama Triple Negativo (TNBC)

Repositorio con los scripts y matrices de trabajo utilizados en el Trabajo Final de Máster (UOC):
**“Análisis integrativo multiómico en cáncer de mama triple negativo para la identificación de perfiles moleculares asociados al pronóstico clínico”**.

---

## 🎯 Objetivo

Documentar y disponibilizar el pipeline bioinformático para el análisis multiómico en TNBC, incluyendo:

- Transcriptómica (RNA-seq)
- Metilación de ADN
- Copy Number Variation (CNV)
- Single Nucleotide Variations (SNV)
- Integración multiómica (MOFA2) y generación de resultados/figuras

---

## 📂 Contenido del repositorio

### Scripts (R)
Pipeline en 6 etapas:

1. `code_1_clinical_data.R`  
   Descarga/procesamiento de datos clínicos y filtrado de cohorte TNBC.
2. `code_2_transcriptomic_layer.R`  
   Preprocesamiento/análisis de la capa transcriptómica.
3. `code_3_methylation_layer.R`  
   Preprocesamiento/análisis de metilación de ADN.
4. `code_4_snv_layer.R`  
   Preprocesamiento/análisis de SNV (mutaciones somáticas).
5. `code_5_cnv_layer.R`  
   Preprocesamiento/análisis de CNV.
6. `code_6_integration.R`  
   Integración multiómica (MOFA2) + outputs.

### Matriz de con seleccion final de pacientes TNBC (xlsx)
Matriz final (n = 60) con los identificadores de pacientes con fenotipo TNBC que cuentan con información disponible en todas las capas ómicas analizadas.

- `TNBC_Patients_Final.xlsx` (cohorte final / IDs / selección final)

### Matrices de integración (xlsx)
Matrices procesadas y exportadas para la integración multiómica, correspondientes a cada una de las capas de datos:

- `TNBC_For.Integration_Clinical.xlsx`
- `TNBC_For.Integration_Transcriptomics.xlsx`
- `TNBC_For.Integration_Methylation.xlsx`
- `TNBC_For.Integration_CNV.xlsx`
- `TNBC_For.Integration_SNV.xlsx`
## Copyright and Code Reuse

Copyright © 2026 Zanie Stella Pankow Farina. All Rights Reserved.

This repository is publicly available for academic transparency, review, and reproducibility purposes only.

No permission is granted to copy, reuse, modify, redistribute, sublicense, commercialize, or incorporate the source code, in whole or in substantial part, into another work without prior written permission from the copyright holder.

Any authorized reuse of the source code must retain attribution to Zanie Stella Pankow Farina and cite the associated work.

For permission to reuse the source code, please contact the copyright holder.

### Citation

Pankow, Z. (2026). *Integrative multi-omics analysis of triple-negative breast cancer to identify molecular profiles associated with clinical prognosis*. Zenodo. DOI: 10.5281/zenodo.19791106
---


