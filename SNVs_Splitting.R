#### Function: split MNVs (multiple-nucleotide variants) into different SNVs (single-nucleotide variants). Output a new dataframe containing
#              the SNVs.

# variants.table: dataframe containing the information from the filtered SNVs in our standard format.
# context: if TRUE, the sequence context for each splitted SNV is retrieved from the reference (-10bp, +11bp). It is much slower.
# genome.ref: DNAStringSet object (from BioStrings) containing the Tasmanian devil reference genome 7.1. Needed only if the sequence context needs to be retrieved.

#              library('Biostrings');
#              genome <- readDNAStringSet(filepath='/Users/dmh/Desktop/utilities/genomes/7.1/genome.fa'); 


split.mnvs <- function(variants.table, context=FALSE, genome.ref=NA){
  
  ## Extract those rows which contain MNVs.
  
  select <- nchar(as.character(variants.table[,3])) > 1;
  mnvs <- variants.table[select,];
  
  ## Split the MNVs.
  
  print('Splitting the MNVs ...');
  
  # Initialize empty dataframe to contain the splitted MNVs.
  
  split.mnvs <- as.data.frame(setNames(replicate(ncol(mnvs),numeric(0), simplify = F), colnames(mnvs)));
  
  # Obtain splitted MNVs without the sequence context.
  
  split.mnvs <- apply(mnvs, 1, function(x){
    
    # Compare reference and altered sequences and extract the positions that are different.
    
    ref.seq <- as.character(unlist(x[3]));
    alt.seq <- as.character(unlist(x[4]));
    changed.positions <- which(unlist(strsplit(ref.seq, '')) != unlist(strsplit(alt.seq, '')));
    
    # Extract information needed for the splitted MNVs (columns 1-6, 8).
    
    contig <- unlist(strsplit(as.character(unlist(x[1])), ':'))[1];
    con.positions <- as.numeric(unlist(strsplit(as.character(unlist(x[1])), ':'))[2]) + (changed.positions-1);
    col.1 <- paste0(contig, ':', con.positions); # Column 1
    
    col.2 <- as.numeric(unlist(x[2])) + (changed.positions-1); # Column 2
    
    col.3 <- unlist(strsplit(ref.seq, ''))[changed.positions]; # Column 3
    
    col.4 <- unlist(strsplit(alt.seq, ''))[changed.positions]; # Column 4
    
    col.5 <- rep(as.character(unlist(x[5])), length(changed.positions)); # Column 5
    
    col.6 <- rep(as.character(unlist(x[6])), length(changed.positions)); # Column 6
    
    col.7 <- rep('', length(changed.positions)); # Column 7, only initialize
    
    col.8 <- rep(as.character(unlist(x[8])), length(changed.positions)); # Column 8
    
    return(cbind(col.1,col.2,col.3,col.4,col.5,col.6,col.7,col.8));
    
  })
  
  split.mnvs.final <- do.call(rbind, split.mnvs);
  
  
  # Add sequence context if required.
  
  if(context){
    
    print('Adding sequence context to the splitted MNVs ...');
  
    all.contigs <- sapply(strsplit(split.mnvs.final[,1], ':'), '[[', 1);
    all.contigs.pos <- sapply(strsplit(split.mnvs.final[,1], ':'), '[[', 2);
    info.context <- cbind(all.contigs, all.contigs.pos);
  
    seq.contexts <- unlist(apply(info.context, 1, function(x){
    
      genome.ref[x[1]][[1]][(as.numeric(x[2])-10):(as.numeric(x[2])+10)]
    
    }));
    
    seq.contexts.char <- sapply(seq.contexts, as.character);
    split.mnvs.final[,7] <- seq.contexts.char;
    
  }
  
  # Finish formatting of splitted MNVs.
  
  split.mnvs.final <- as.data.frame(split.mnvs.final);
  names(split.mnvs.final) <- names(variants.table);
  
  
  ## Merge together the original SNVs and the new splitted SNVs. Avoid duplicates and sort by chromosomal position.
  
  print('Creating final dataframe ...');
  
  # Extract original SNVs.
  
  original.snvs <- variants.table[nchar(as.character(variants.table[,3])) == 1,];
  
  # Merge all the SNVs (original, splitted).
  
  all.snvs <- rbind(original.snvs, split.mnvs.final);
  
  # Sort and remove possible duplicates.
  
  all.snvs.final <- unique(all.snvs[order(as.numeric(all.snvs[,2])),]);
  
  
  ## Return the final dataframe.
  
  return(all.snvs.final);
  
}
