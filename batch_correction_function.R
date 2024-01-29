# ----------------------------------------------------------------------------------
# Title: Batch Correction function for HER2 CN RNA Seq Data
# Author: Hina Dalal
# Date of Analysis: 2023-09-04
# Replace with `Sys.Date()` for the current date
# Description: Gene expression values will be normalized to reference batch (refbatch) by mean centering of log2 transformed TPM GEX matrix
# to derive gene specific mean differences between library protocols. The calculated differences (scaling.factors) were used as correction factor 
# for converting gene expression between protocols.
# Version: 1.0 
# Project Title: Quantitative digital PCR measurement of ERBB2 copy number is predictive of outcome in early breast cancer 
# patients treated with adjuvant trastuzumab
# ----------------------------------------------------------------------------------



# Actual R code starts below



######################################################################
########### Batch Correction function for HER2 CN RNA Seq Data #######
######################################################################

# Make sure that genes (features) are in rows and samples in col in gex_mat: GEX (Gene Expression Matrix)
# batches - sample annotations containing information on batches
# Make sure batches is a factor with atleast two or more levels and has same length as ncol(gex_mat)

adjustmean.ref <- function(gex_mat,batches,refbatch) {
    
    # some checks first
    if (!is.factor(batches) || length(levels(batches)) < 2) {
      stop("batches is not a factor or has fewer than two levels")
    }
    
    if (any(is.na(batches))) {
      stop("batches contains NA values")
    }
    
    if (length(batches) != ncol(gex_mat)) {
      stop("batches has not the same length as ncol(gex_mat)")
    }
    
    if (any(table(batches) < 1.5)) {
      stop("a level of batches has one or zero counts")
    }
    
    if (length(refbatch) != 1) {
      stop("refbatch does not have length 1")
    }
    
    if (!is.character(refbatch)) {
      stop("refbatch is not a character")
    }
    
    if (!refbatch %in% levels(batches)) {
      stop("refbatch is not a level of batches")
    }
    
    ## start batch correction
    
    mat_after <- gex_mat
    mat_after[] <- NA
    adjustedvalues <- matrix(ncol=length(levels(batches)),nrow=nrow(gex_mat),dimnames=list(rownames(gex_mat),levels(batches)))
    refindex <- which(batches==refbatch)
    mat_ref <- gex_mat[,refindex]
    rowMeanref <- apply(mat_ref, 1, mean, na.rm = T)
    
    
    for (i in 1:length(levels(batches))){
      index <- which(batches==levels(batches)[i])  
      mat1 <- gex_mat[,index]
      
      mat1adj <- mat1
      mat1adj[] <- NA
      for (j in 1:nrow(mat1)){
        adjusted <- rowMeanref[j] - mean(mat1[j,],na.rm=T)
        adjustedvalues[j,i] <- adjusted
        mat1adj[j,] <- mat1[j,]+adjusted
      }
      mat_after[,index] <- mat1adj
    }
    return(list(adjusted.data=mat_after,scaling.factors=adjustedvalues))
  }

### Example Usage

# erbb2_adj <-adjustmean.ref(log2(TPM_mat+0.1),clinical_df$Protocol,  "dUTP") 


