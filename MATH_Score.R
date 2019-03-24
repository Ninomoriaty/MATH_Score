# testing
# library(maftools)

# maf_file = "/home/ninomoriaty/R_Project/311252_snv_indel.imputed.maf"
######## MATH Score ##########
MATH_score <- function(maf_file, tsb = "OFA", minvaf = 0, maxvaf = 1){
  # read .maf file
  maf_input <- read.table(maf_file, quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # get vaf-related infomation
  dat.hugo_symbol <- maf_input[,1]
  dat.vaf <- maf_input[,ncol(maf_input)-3]
  dat.tsb <- data.frame(maf_input[,ncol(maf_input)])
  vaf_input_mt <- data.frame(dat.hugo_symbol, dat.vaf, dat.tsb)
  colnames(vaf_input_mt) <- c("Hugo_Symbol", "VAF", "Tumor_Sample_Barcode")
  # get all sample names
  samples_math <- data.frame()
  samples <- data.frame(maf_input[,ncol(maf_input)])
  tsb_ls <- as.data.frame(as.data.frame(table(samples))["samples"][which(as.data.frame(table(samples))["samples"]$samples != ""),])
  
  # MATH Result
  if (tsb == "OFA"){
    # list all samples' MATH scores
    for (counter_mt in 1:length(tsb_ls[,1])){
      for (sample_name_mt in tsb_ls){
        VAF_column <- data_clean(vaf_input_mt, as.character(sample_name_mt)[counter_mt], minvaf, maxvaf)
        sample_math <- data.frame(as.character(sample_name_mt)[counter_mt], math_cal(VAF_column))
        samples_math <- rbind(samples_math, sample_math)
      }
    }
    colnames(samples_math) <- c("Tumor_Sample_Barcode", "MATH_score")
    samples_math
  } else{
    # calculate specific sample's MATH score
    VAF_column <- data_clean(vaf_input_mt, tsb, minvaf, maxvaf)
    math_cal(VAF_column)
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
  MAD = median(abs(VAF_column - median(VAF_column)))
  MATH = 100 * MAD * 1.4826/median(VAF_column)
  MATH
}


######### Compare with maftools ###############
# laml = read.maf(maf = maf_file)
# tcga.ab.2972.het = inferHeterogeneity(maf = laml, tsb = '311252-V', vafCol = 'VAF', useSyn = TRUE)
# VAF_column2 <- tcga.ab.2972.het$"clusterData"$t_vaf
# plotClusters(clusters = tcga.ab.2972.het)








