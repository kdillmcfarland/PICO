"
Nonparametric Tests for Repeated Measures Data in Factorial Designs across multiple taxa

x = Otu table
taxa.list = OTUs, genera, etc of interest in the form c('','',''...)
output_name = name of ouput file
"

nparLD.pval.fxn = function(x, taxa.list, output_name)
{require(nparLD)
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
      # Spread and then gather to force NAs for all missing measurements for individual subjects
      # Required to make `f1.ld.f1` work since it assume length = subjects x time points
      spread(key=Time_Point, value=taxa) %>% 
      gather(key=Time_Point, value=Taxon, 3:5)
    
    # Run F1-LD-F1 design test
    model_1 = 
      f1.ld.f1(taxa.temp$Taxon, time=taxa.temp$Time_Point, group=taxa.temp$Treatment, subject=taxa.temp$Subject_ID, plot.RTE=FALSE, description=FALSE, order.warning=FALSE)
    
    {
      # pull out p-vals
      # Wald-type test statistic
      wald.p = as.matrix(model_1$Wald.test)
      # ANOVA-type test statistic with Box approximation
      aov.p = as.matrix(model_1$ANOVA.test)
      # modified ANOVA-type test statistic with Box approximation
      aov.mod.p = as.matrix(model_1$ANOVA.test.mod.Box)
      # Wald-type test statistic for simple time effect
      wald.time.p = as.matrix(model_1$Wald.test.time)
      # ANOVA-type test statistic for simple time effect
      aov.time.p = as.matrix(model_1$ANOVA.test.time)
      
      # Bind pvals into 1 table
      p.vals = data.frame(
        pval = c(wald.p[,"p-value"], 
                 aov.p[,"p-value"], 
                 aov.mod.p[,"p-value"], 
                 wald.time.p[,"p-value"], 
                 aov.time.p[,"p-value"]),
        test = c(rep("wald", 3), 
                 rep("aov", 3), 
                 rep("aov.mod", 1), 
                 rep("wald.time", 2), 
                 rep("aov.time", 2)),
        variable = c(rownames(wald.p),
                     rownames(aov.p),
                     rownames(aov.mod.p),
                     rownames(wald.time.p),
                     rownames(aov.time.p)),
        taxon = rep(taxa, 11))
      
      p.vals=as.matrix(p.vals)
      
      write.table(p.vals, file=paste(output_name, ".nparLD.csv", sep=""), append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)
      }}
  
  y=read.table(paste(output_name, ".nparLD.csv", sep=""), header=FALSE, sep=",", fill = TRUE)
  file.remove(paste(output_name, ".nparLD.csv", sep = ""))
  colnames(y) = c("pval", "test", "variable", "taxon")
  write.csv(y,file=paste(output_name,'.nparLD.csv', sep=''), row.names=FALSE)
}