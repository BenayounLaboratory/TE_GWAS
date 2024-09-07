Run_ORA <- function(output.dir, output_subfolder, universe.list, my.genelist, my.gs.collection, gs.label, condition_label, output.object) {
  
  
  # Define the output directory
  dir.output <- paste(output.dir, output_subfolder, sep = '')
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT
  # enricher uses a weird universe, which is defined as the intersection of the supplied universe list with gene set collection genes. This is problematic when the gene sets don't cover the entire genome
  # Here, add the supplied universe list to the gene set collection. The intersection of the supplied universe list with the gene sets should now be the supplied universe list (which is the background we want)
  # Adjust the max GSSize so it does not include the universe geneset
  
      # Create a df to hold the universe geneset
      universe.geneset <- data.frame(gs_name = 'Universe', gene = universe.list)
      
      # Update the colnames to match the geneset collection
      colnames(universe.geneset) <- colnames(my.gs.collection)
      
      # Add the universe geneset to the geneset collection of interest
      my.gs.collection <- rbind(my.gs.collection, universe.geneset)
  
  # Run ORA
  my.ORA <- enricher(gene = my.genelist,
                     universe      = universe.list,
                     TERM2GENE     = my.gs.collection,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1,
                     minGSSize     = 15,
                     maxGSSize     = 10000)
  
  # if no genes mapped to the gene list, end the function
  if (is.null(my.ORA)) {
    return()
  }
  
  # Save *ALL* results, regardless of significance
  write.table(my.ORA@result, file = paste(dir.output, 'ClusterProfiler_ORA_Table_', gs.label, '_', condition_label, '_ALL', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
  
  # Make an object to hold results filtered on significance
  my.ORA.sig <- my.ORA
  
  # Filter results on FDR < 0.05
  my.ORA.sig@result <- my.ORA.sig@result[which(my.ORA.sig@result$p.adjust < 0.05 & my.ORA.sig@result$qvalue < 0.20), ]
  
  # Define the number of significant results
  my.sig.gs.num <- nrow(my.ORA.sig@result)
  
  # if there are significant results, save files
  if (my.sig.gs.num > 0) {

      # Define the number of terms to plot. Do a maximum of 10 *significant* terms
      if (my.sig.gs.num > 10) {
        terms_to_plot <- 10
      } else {
        terms_to_plot <- my.sig.gs.num
      }
      
      # Define the maximum gene ratio
      generatio_limit <- ceiling(max(as.numeric(my.ORA.sig@result[1:terms_to_plot, 'GeneRatio'])))
      
      # write results to file
      write.table(my.ORA.sig@result, file = paste(dir.output, 'ClusterProfiler_ORA_Table_', gs.label, '_', condition_label, '_FDR5_', my.sig.gs.num, '_Significant' ,'.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
  
      # Make and save dotplot with top genesets
      pdf(paste(dir.output, 'ClusterProfiler_ORA_Dotplot_', gs.label, '_', condition_label, '_FDR5', ".pdf", sep=""), width = 9, height = 8)
  
        print(dotplot(my.ORA.sig, x = TRUE, showCategory = terms_to_plot, size = '-log10(p.adjust)', color = NULL, title = paste(condition_label, " ORA ", gs.label, sep = ''), font.size = 12) + scale_size(range = c(4, 12)) + aes(shape = I(16)) + aes(color = `GeneRatio`) + set_enrichplot_color(type = "color", name = "GeneRatio", colors = c("white","red"), limits=c(-0.01, generatio_limit))) 
  
      dev.off()

  }
  
  # No output needed
  return()
  
  
} # END FUNCTION
