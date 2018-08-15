"
Dunn's Kruskal-Wallis Multiple Comparisons across multiple taxa

x = Otu table
taxa.list = OTUs, genera, etc of interest in the form c('','',''...)
output_name = name of ouput file
"

dunn.pval.fxn = function(x, taxa.list, output_name)
{require(FSA)
 require(tidyverse)
  
  for(taxa in taxa.list)
  {
    # Subset metadata to needed variables
    meta.temp = meta %>% 
      rownames_to_column() %>% 
      dplyr::select(rowname, Subject_ID, Time_Point, Treatment)
    
    # Pull out taxon of interest and format for statistics test
    taxa.temp = x %>% 
      rownames_to_column() %>% 
      # Pull out taxon of interest
      dplyr::select(rowname, taxa) %>% 
      # Join taxon counts with subset metadata
      full_join(meta.temp, by="rowname") %>% 
      dplyr::select(-rowname) %>% 
      mutate(group=as.factor(paste(Time_Point, Treatment, sep="_"))) %>% 
      rename(Taxon = taxa)
    
    # Run F1-LD-F1 design test
    model_1 = dunnTest(taxa.temp$Taxon ~ taxa.temp$group, method="bh")
    
    {
      # pull out p-vals
      p.vals = data.frame(
        fdr_pval = model_1$res$P.adj,
        label = model_1$res$Comparison,
        taxon = rep(taxa, length(model_1$res$P.adj)))
      
      p.vals=as.matrix(p.vals)
      
      write.table(p.vals, file=paste(output_name, ".dunn.fdr.csv", sep=""), append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)
    }}
  
  y=read.table(paste(output_name, ".dunn.fdr.csv", sep=""), header=FALSE, sep=",", fill = TRUE)
  file.remove(paste(output_name, ".dunn.fdr.csv", sep = ""))
  colnames(y) = c("fdr_pval", "label", "taxon")
  write.csv(y,file=paste(output_name,'.dunn.fdr.csv', sep=''), row.names=FALSE)
}