# testing
#library(maftools)

#maf_file = "/home/ninomoriaty/R_Project/311252_snv_indel.imputed.maf"
######## MATH Score ##########
MATH_score <- function(maf_file){
  # read.maf
  maf_input <- read.table(maf_file, header = TRUE, fill = TRUE)
#  samples <- data.frame(maf_input[,ncol(maf_input)])
#  tsb_ls <- as.data.frame(as.data.frame(table(samples))["samples"][which(as.data.frame(table(samples))["samples"]$samples != ""),])
  
  # Data cleaning
  VAF_column = maf_input[which(maf_input$Tumor_Sample_Barcode == "311252-S"),]$VAF
  VAF_column = VAF_column[which(!is.na(VAF_column))]
  VAF_column = as.numeric(as.character(VAF_column))
  
  # MATH Caculation
  MAD = median(abs(VAF_column - median(VAF_column)))
  MATH = 100 * MAD * 1.4826/median(VAF_column)
  
}


######### Compare with maftools ###############
# laml = read.maf(maf = maf_file)
# tcga.ab.2972.het = inferHeterogeneity(maf = laml, tsb = '311252-S', vafCol = 'VAF')
# VAF_column <- tcga.ab.2972.het$"clusterData"$t_vaf
# plotClusters(clusters = tcga.ab.2972.het)








