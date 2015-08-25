############ SNP_Prepare - DFT2 Mutational Signatures ##############
########################## Version 1.0 #############################

## INPUT AND PREPARE FUNCTIONS
## Last Update - 24/08/2015 ##

library(stringr)

# Set Directories
Home <- c("/Users/ms37/Desktop/Data/Platypus_Full")
setwd(Home)

# Prepare the two raw tumour files
files <- list.files()[grep("#", list.files())]
for (i in c(2,4)){
  
  print(files[i])
  
  ## a. Load Platypus-File, takes ~ 10 minutes
  setwd("/Users/ms37/Desktop/Data/Platypus_Full/")
  platypus.name1 <- list.files()[i]
  tmp <- file(paste0(platypus.name1))
  tmp.full <- read.table(files[i])[,-c(3,9)]
  colnames(tmp.full) <- c("Contig", "Position", "Base after deletion/before insertion", "Pattern", "Quality", "Filter", "Stats", "Tumour", "Normal")    
  
  ## b. Trimming: Remove non-Deletions/Insertions
  
  # Lengths of Strings
  ref.lengths <- nchar(as.character(tmp.full[,3]))
  alt.lengths <- nchar(as.character(tmp.full[,4]))
  length.check <- ref.lengths==alt.lengths
  tmp.trimmed <- tmp.full[which(length.check==TRUE),]
  
  ## c. Split up into the different Chromosomes
  setwd(paste0("/Users/ms37/Desktop/Data/SNS_PLATYPUS/", files[i]))
  for (j in c(1,2,3,4,5,6,"U","x")){
    tmp <- tmp.trimmed[grep(paste0('Chr',j), tmp.trimmed[,'Contig']),]
    
    # Output
    write.csv(tmp, paste0("platy_", strsplit(files[i], "_")[[1]][1], "_chr", j,".csv"), row.names=F)
  }
  rm(tmp, tmp.full, tmp.trimmed, ref.lengths, alt.lengths, j, length.check)
}

# Remove Normal column, extract SC
contig.size <- read.table("/Users/ms37/Desktop/Data/Info-Files/Ref7.1_ContigPositions.txt", header=T)
check <- read.csv("platy_#203T3_chrU.csv")
check.stats <- str_split_fixed(check[,"Stats"],";", 13)[,12]
check.stats <- str_split_fixed(check.stats,"=", 2)[,2]
check[,"Stats"] <- check.stats
check <- check[,-9]
check[,1] <- paste(check[,1], check[,2], sep=":")
check.calls <- str_split_fixed(check[,"Tumour"],":", 6)
check.calls <- paste(check.calls[,6], check.calls[,5], sep=":")
check.calls <- str_split_fixed(check.calls, ":", 2)
check.calls[grep("[,]", check.calls[,1]),1] <- str_split_fixed(check.calls[grep("[,]", check.calls[,1]),1], "[,]", 2)[,1]
check.calls[grep("[,]", check.calls[,2]),2] <- str_split_fixed(check.calls[grep("[,]", check.calls[,2]),2], "[,]", 2)[,1]
check.calls <- paste(check.calls[,1], check.calls[,2],
                     as.integer(check.calls[,1])/as.integer(check.calls[,2]), sep=":")
check[,"Tumour"] <- check.calls
check.contig <- str_split_fixed(check[,1], ":", 2)[,1]
check[which(check.contig!="ChrU_supercontig_000000000"),2] <- 
  as.integer(check[which(check.contig!="ChrU_supercontig_000000000"),2])+
  contig.size[match(check.contig[which(check.contig!="ChrU_supercontig_000000000")],contig.size[,1]),3]
SN.Tumour.ChrU.unfiltered <- check
rm(list=ls()[-6])
save.image("SN.Tumour.ChrU.unfiltered.Rdata")
rm(list=ls())
