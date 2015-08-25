############ SNP_Overlap - DFT2 Mutational Signatures ##############
########################## Version 1.0 #############################

## OVERLAP AND FILTER FUNCTIONS
## Last Update - 25/08/2015 ##

DFT2_SNVs_Overlap <- function(chromosome, threshold, flag){
 
  require(VennDiagram)
  setwd("/Users/ms37/Desktop/Data/SNS_PLATYPUS/Filtersets")
  
  # Load in chromosomal sets
  load(paste0("/Users/ms37/Desktop/Data/SNS_PLATYPUS/#202T2_SimReads_full.vcf/RV.Tumour.Chr", chromosome, ".unfiltered.Rdata"))
  load(paste0("/Users/ms37/Desktop/Data/SNS_PLATYPUS/#203T3_SimReads_full.vcf/SN.Tumour.Chr", chromosome, ".unfiltered.Rdata"))
  
  RV.T <- get(paste0("RV.Tumour.Chr", chromosome, ".unfiltered"))
  SN.T <- get(paste0("SN.Tumour.Chr", chromosome, ".unfiltered"))
                          
  
  # Take positional intersections and unique sets
  conc <- length(intersect(RV.T[,1], SN.T[,1]))/
    length(union(RV.T[,1], SN.T[,1]))
  
  cat("\n Unfiltered Tumour Concordance: ", round(conc,2), "%")  
  
  
  # Build Venn
  full.RV <- length(RV.T[,1]) 
  full.SN <- length(SN.T[,1])
  intersection <- intersect(RV.T[,1], SN.T[,1])
  intersection <- length(intersection)
  
  png(paste0("DFT2_SNVs_chr_", chromosome, "_prefilter.png"), width=400,height=350)
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 5, 0, 0))
  draw.pairwise.venn(area1=full.RV,
                     area2=full.SN,
                     cross.area=intersection,
                     category = c(paste0("RV - 202T2"),
                                  paste0("SN - 203T3")),
                     fontfamily = rep("sans", 3), cat.fontfamily = "sans",
                     fill = c("blue", "orange"), euler.d=TRUE, scaled=TRUE, 
                     ind = T, cex=1.0, cat.cex=0.7, cat.dist=0.03)
  dev.off()
  
  
  # Build Lists
  
  shared <- intersect(RV.T[,1], SN.T[,1])
  shared.RV <- RV.T[match(shared, RV.T[,1]),]
  shared.SN <- SN.T[match(shared, SN.T[,1]),]
  shared.T <- list(shared.RV,shared.SN)
  names(shared.T) <- c("RV", "SN")
  
  unique.RV <- RV.T[-match(shared, RV.T[,1]),]
  unique.SN <- SN.T[-match(shared, SN.T[,1]),]
  unique.T <- list(unique.RV, unique.SN)
  names(unique.T) <- c("RV", "SN")
  
  
  # Apply quality filters: QUAL field and flags

  shared.T$'RV' <- shared.T$'RV'[shared.T$'RV'[,5]>threshold,]
  shared.T$'SN' <- shared.T$'SN'[shared.T$'SN'[,5]>threshold,]
  
  unique.T$'RV' <- unique.T$'RV'[unique.T$'RV'[,5]>threshold,]
  unique.T$'SN' <- unique.T$'SN'[unique.T$'SN'[,5]>threshold,]
  
  if (flag=="Y"){
    
    shared.T$'RV' <- shared.T$'RV'[shared.T$'RV'[,"Filter"]=="PASS",]
    shared.T$'SN' <- shared.T$'SN'[shared.T$'SN'[,"Filter"]=="PASS",]
   
    unique.T$'RV' <- unique.T$'RV'[unique.T$'RV'[,"Filter"]=="PASS",]
    unique.T$'SN' <- unique.T$'SN'[unique.T$'SN'[,"Filter"]=="PASS",]   
    
  }
  
  # Build Venn
  full.RV2 <- length(unique.T$'RV'[,1])+length(shared.T$'RV'[,1])
  full.SN2 <- length(unique.T$'SN'[,1])+length(shared.T$'SN'[,1])
  intersection2 <- length(shared.T$'SN'[,1])
  
  png(paste0("DFT2_SNVs_chr_", chromosome, "_postfilter.png"), width=400,height=350)
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 5, 0, 0))
  draw.pairwise.venn(area1=full.RV2,
                     area2=full.SN2,
                     cross.area=intersection2,
                     category = c(paste0("RV - 202T2"),
                                  paste0("SN - 203T3")),
                     fontfamily = rep("sans", 3), cat.fontfamily = "sans",
                     fill = c("blue", "orange"), euler.d=TRUE, scaled=TRUE, 
                     ind = T, cex=1.0, cat.cex=0.7, cat.dist=0.03)
  dev.off()
  
  conc2 <- intersection2/(full.RV2+length(unique.T$'SN'[,1]))
  
  cat("\n Filtered Tumour Concordance: ", round(conc2,2), "%")  
  
  
  # Output: Filtered Unique and Intersection Lists
  
  sets <- list(shared.T$'RV', shared.T$'SN', unique.T$'RV', unique.T$'SN')
  names(sets) <- c("shared.RV", "shared.SN", "unique.RV", "unique.SN")
  return(sets)
}

for (i in c(1,2,3,4,5,6,"x","U")){
  assign(paste0("Tumours.chr",i), DFT2_SNVs_Overlap(chromosome = i, threshold = 150, flag = "Y"))
}