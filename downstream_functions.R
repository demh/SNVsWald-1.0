##### Downstream analysis of SNVs #######

#### Function: create a mutational signature for the SNVs identified. STILL NEED TO CODE FOR CONTEXT

# variants.table: dataframe containing the information from the filtered SNVs in our standard format.
# context: if TRUE, sequence context is taken into account to construct the signature.
# output.path: absolute path to output the plot.
# title.plot: main title in the mutational signature plot.

mut.signature <- function(variants.table, context=FALSE, output.path, title.plot=''){
  
  ## Pyrimidine-context substitutions.
  
  equivalent.mut <- matrix(c('C>A', 'G>T', 
                             'C>T', 'G>A', 
                             'C>G', 'G>C', 
                             'T>A', 'A>T', 
                             'T>C', 'A>G', 
                             'T>G', 'A>C'),
                           nrow=6, ncol=2, byrow=TRUE);
  
  ## Take into account sequence context.
  
  if(context){
    
    
  }else{ ## Do not take into account sequence context
    
    # Extract the variants from the input table and save them in consensus format (e.g. C>T).
    
    snvs <- paste0(variants.table[,3], '>', variants.table[,4]);
    
    # Convert to pyrimidine context.
    
    for(mut in 1:nrow(equivalent.mut)){
      
      snvs <- gsub(equivalent.mut[mut,2], equivalent.mut[mut,1], snvs);
    }
    
    # Absolute counts of each type of SNVs.
    
    consensus.mut <- equivalent.mut[,1];
    counts <- sapply(consensus.mut, function(x){length(grep(x, snvs))});
    
    # Plot the results.
    
    pdf(output.path, height=10, width=10);
    
    barplot(counts, col=c('blue', 'red', 'black', 'grey', 'green', 'magenta'),
            main=title.plot, cex.main = 0.9);
    
    dev.off(); 
    
  }
}


#### Function: create a BAF plot for the SNVs in one chromosome. The Y-axis represents the BAF values and the X-axis the chromosomal positions.

# variants.table: dataframe containing the information from the filtered SNVs in our standard format.
# colour: if TRUE, colour the SNVs in the plot according to the type of substitution (pyrimidine consensus).
# split: if TRUE, one plot per each type of mutation is created.
# output.path: absolute path to output the plot.
# title.plot: main title in the mutational signature plot.

baf.plot <- function(variants.table, colour=TRUE, split=FALSE, output.path, title.plot=''){
  
  ## Extract the information which is needed (chromosomal position, BAF, type of substitution).
  
  baf.info <- sapply(strsplit(variants.table[,8], ':'), '[[', 3);
  useful.info <- cbind(variants.table[,2], baf.info, paste0(variants.table[,3], '>', variants.table[,4]));
  
  ## Convert 'type of substitution' to pyrimidine context.
  
  # Pyrimidine-context substitutions.
  
  equivalent.mut <- matrix(c('C>A', 'G>T', 
                             'C>T', 'G>A', 
                             'C>G', 'G>C', 
                             'T>A', 'A>T', 
                             'T>C', 'A>G', 
                             'T>G', 'A>C'),
                           nrow=6, ncol=2, byrow=TRUE);
  
  # Conversion.
  
  for(mut in 1:nrow(equivalent.mut)){
    
    useful.info[,3] <- gsub(equivalent.mut[mut,2], equivalent.mut[mut,1], useful.info[,3]);
  }
  
  
  ## Add a fourth column to the useful.info with the conversion from 'type of substitution' to colour for the plotting.
  
  colours.info <- useful.info[,3];
  colours <- c('blue', 'red', 'black', 'grey', 'green', 'magenta');
  
  for(col in 1:length(colours)){
    
    colours.info <- gsub(equivalent.mut[col,1], colours[col], colours.info);
  }
  
  useful.info <- cbind(useful.info, colours.info);
  
  
  ## Plot the results.
  
  if(colour){
  
    if(split){
      
      pdf(output.path, height=20, width=24);
      
      par(mfrow=c(3,2));
      
      for(mut in 1:length(equivalent.mut[,1])){
        
        select <- useful.info[,3] == equivalent.mut[mut,1];
        
        plot(as.numeric(useful.info[select,1]), as.numeric(useful.info[select,2]), col=useful.info[select,4], main=paste0(title.plot,'_', equivalent.mut[mut,1]), 
             pch='.', xlab='Chromosomal position', ylab='BAF', cex.axis=0.8, ylim=c(0,1));
        legend('topright', legend=equivalent.mut[mut,1], col=colours[mut], pch=16, cex=0.8, ncol=2, bg='white');
        
      }
      
      dev.off();
      
    }else{
    
      pdf(output.path, height=6, width=12);
    
      plot(as.numeric(useful.info[,1]), as.numeric(useful.info[,2]), col=useful.info[,4], main=title.plot, pch='.',
        xlab='Chromosomal position', ylab='BAF', cex.axis=0.8, ylim=c(0,1));
      legend('topright', legend=equivalent.mut[,1], col=colours, pch=16, cex=0.8, ncol=2, bg='white');
    
      dev.off();
    }
  
  }else{
    
    pdf(output.path, height=6, width=12);
    
    plot(as.numeric(useful.info[,1]), as.numeric(useful.info[,2]), main=title.plot, pch='.',
         xlab='Chromosomal position', ylab='BAF', cex.axis=0.8, ylim=c(0,1));
    
    dev.off(); 
  }  
}
