# HER2_CNA_ddPCR
RNA Seq data analysis code for batch correction for ddPCR HER2 copy number amplification project
## Overview

This repository contains an R script for batch correction in gene expression data, particularly for HER2 CN RNA Sequencing data. The function normalizes gene expression values to a reference batch (refbatch) by mean centering of log2 transformed TPM GEX matrix to derive gene-specific mean differences between library protocols. The calculated differences (scaling factors) are used as correction factors for converting gene expression between protocols.

### Project Title
Quantitative digital PCR measurement of ERBB2 copy number is predictive of outcome in early breast cancer patients treated with adjuvant trastuzumab.

## Author
Hina Dalal

## Date of Analysis
2024-01-29  
(Replace with `Sys.Date()` for the current date)

## Version
1.0 

## Description
The provided R script `adjustmean.ref` function allows for effective batch correction in RNA-Seq data analysis. It is particularly tailored for datasets where gene expression needs to be standardized across different batch protocols. By mean centering and applying calculated scaling factors, the function ensures comparability and consistency of gene expression measurements across different batches.

## Usage

To use the function, ensure that your gene expression matrix (`gex_mat`) has genes (features) in rows and samples in columns. The `batches` variable should be a factor with at least two levels and the same length as the number of columns in `gex_mat`.

### Example
```R
# erbb2_adj <- adjustmean.ref(log2(TPM_mat+0.1), clinical_df$Protocol, "dUTP")
