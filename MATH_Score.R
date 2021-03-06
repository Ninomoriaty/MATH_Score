#' A Multifunctional MATH Calculator
#' @description read maf file as data.frame, calculate MATH score and present results in different ways determined by parameters.
#' This function requires VAF for MATH score calculation and Tumor_sample_Barcode for selection.
#' VAF can be on the scale 0-1 and you could filter VAF with minvaf and maxvaf.
#' 
#' @param maf_file specify a maf document/directory as the input of the function
#' @param tsb specify single/all sample names (Tumor_Sample_Barcodes)
#' @param minvaf filter low frequency variants caused by sequencing error. Default 0. (on the scale of 0 to 1)
#' @param maxvaf filter high frequency variants due to copy number alterations or impure tumor. Default 1. (on the scale of 0 to 1)
#' @return MATH score list for all/selected samples or a single MATH score value

#' @examples
#' \dontrun{
#' maf_file = "/home/ninomoriaty/R_Project/311252_snv_indel.imputed.maf"
#' MATH_score(maf_file)
#' MATH_score(maf_file, tsb = c("tsb1"))
#' MATH_score(maf_file, tsb = c("tsb1", "tsb2", "tsb3"))
#'}


######## MATH Score ##########
MATH_score <- function(maf_file, tsb = c("OFA"), minvaf = 0, maxvaf = 1){
  # read .maf file
  maf_input <- read.table(maf_file, quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # get vaf-related infomation
  dat.hugo_symbol <- maf_input[,1]
  dat.vaf <- maf_input[,ncol(maf_input)-3]
  dat.tsb <- data.frame(maf_input[,ncol(maf_input)])
  vaf_input_mt <- data.frame(dat.hugo_symbol, dat.vaf, dat.tsb)
  colnames(vaf_input_mt) <- c("Hugo_Symbol", "VAF", "Tumor_Sample_Barcode")
  # get all sample names
  samples <- data.frame(maf_input[,ncol(maf_input)])
  tsb_ls <- as.data.frame(as.data.frame(table(samples))["samples"][which(as.data.frame(table(samples))["samples"]$samples != ""),])

  # MATH main function
   # MATH results for one/all sample
    if (any(tsb == c("OFA"))){
      # list all samples' MATH scores
      math_ofa <- math_msp(vaf_input_mt, tsb_ls, minvaf, maxvaf)
      print(math_ofa)
      math_all <- math_patient(vaf_input_mt, minvaf, maxvaf)
      cat("ITH MATH score:", as.character(math_all))
    } else{
      # calculate specific samples' MATH score
      tsb_ls <- data.frame(tsb)
      math_sp <- math_msp(vaf_input_mt, tsb_ls, minvaf, maxvaf)
      print(math_sp)
      math_all <- math_patient(vaf_input_mt, minvaf, maxvaf)
      cat("ITH MATH score:", as.character(math_all))
    }
}


# Data cleaning
data_clean <- function(vaf_input_mt, tsb, minvaf, maxvaf){
  VAF_column = vaf_input_mt[which(vaf_input_mt$Tumor_Sample_Barcode == tsb),]$VAF
  VAF_column = VAF_column[which(!is.na(VAF_column))][which(VAF_column > minvaf & VAF_column < maxvaf)]
  VAF_column = as.numeric(as.character(VAF_column))[which(!is.na(VAF_column))]
  VAF_column
}

# MATH Caculation
math_cal <- function(VAF_column){
  MAD_fac = 1.4826*median(abs(VAF_column - median(VAF_column)))
  MATH = 100 * MAD_fac / median(VAF_column)
  MATH
}

# MATH multi-sample process
math_msp <- function(vaf_input_mt, tsb_ls, minvaf, maxvaf){
  samples_math <- data.frame()
  for (counter_mt in 1:length(tsb_ls[,1])){
    for (sample_name_mt in tsb_ls){
      VAF_column <- data_clean(vaf_input_mt, as.character(sample_name_mt)[counter_mt], minvaf, maxvaf)
      sample_math <- data.frame(as.character(sample_name_mt)[counter_mt], math_cal(VAF_column))
      samples_math <- rbind(samples_math, sample_math)
    }
  }
  colnames(samples_math) <- c("Tumor_Sample_Barcode", "MATH_score")
  samples_math
}

# MATH patient calcualtion
math_patient <- function(vaf_input_mt, minvaf, maxvaf){
  VAF_column = vaf_input_mt$VAF
  VAF_column = VAF_column[which(!is.na(VAF_column))][which(VAF_column > minvaf & VAF_column < maxvaf)]
  VAF_column = as.numeric(as.character(VAF_column))[which(!is.na(VAF_column))]
  math_cal(VAF_column)
}

######### Compare with maftools ###############
# testing
# library(maftools)
# laml = read.maf(maf = maf_file)
# tcga.ab.2972.het = inferHeterogeneity(maf = laml, tsb = '311252-V', vafCol = 'VAF', useSyn = TRUE)
# VAF_column2 <- tcga.ab.2972.het$"clusterData"$t_vaf
# plotClusters(clusters = tcga.ab.2972.het)








