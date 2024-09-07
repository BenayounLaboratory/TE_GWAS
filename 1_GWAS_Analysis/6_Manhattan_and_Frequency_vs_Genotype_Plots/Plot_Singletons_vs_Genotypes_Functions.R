plot_SV_genotype_vs_Singleton_Freq <- function(snp_feature_df, sample_info, singleton_column_name, snp.mapping.info, dir.output, output_name, AFR_binary_genotypes, AMR_binary_genotypes, EUR_binary_genotypes, EAS_binary_genotypes, SAS_binary_genotypes, legend_position){
  
    # GENERATE 3X3 PLOTS FOR 9 SVs VS INSERTION FREQUENCY RELATIONSHIPS
  
    # SNP_FEATURE_DF: MAKE SURE ITS ONLY 3 COLUMNS, 1 IS THE SV ID, 2 IS THE ODDS RATIO, and 3 is the FDR
    # SAMPLE_INFO: Contains the number of global singletons in the column specificied by singleton_column_name
    # binary_genotype objects are from BEDMATRIX

  
  
  
    # DEFINE GENERAL PARAMETERS
  
    # Update the colnames for snp_feature_df
    colnames(snp_feature_df) <- c('SNP', 'OR.R.', 'FDR')

    # Define vector of ethnic groups to loop over
    ethnic_groups <- c('AFR', 'AMR', 'EUR', 'EAS', 'SAS')
    
    # Define colors for points
    ancestral.group.colors <- c(rep('#648FFF', 3),
                                rep('#785EF0', 3),
                                rep('#FE6100', 3),
                                rep('#DC267F', 3),
                                rep('#FFB000', 3)
                                )
    
    
  
    
    
    # START THE PLOT OBJECT
    pdf(paste(dir.output, output_name, ".pdf", sep=""), width = 10, height = 10)
    par(mfrow = c(3,3), pty ="m") # 3 rows by 3 columns
    
        # Loop over the entries you would like to plot
        for (ith_relationship in 1:nrow(snp_feature_df)) {
          
          
            # DEFINE VARIABLES FOR THE ITH SV
          
            # Define the ith SV to plot
            SV_to_plot <- snp_feature_df[ith_relationship, 'SNP']
            
                # Define the allele combinations for that SV
                REF_allele <- snp.mapping.info[SV_to_plot, 'REF']
                ALT_allele <- snp.mapping.info[SV_to_plot, 'ALT']
                allele_combo_1 <- paste(REF_allele, '/', REF_allele, sep = '')
                allele_combo_2 <- paste(REF_allele, '/', ALT_allele, sep = '')
                allele_combo_3 <- paste(ALT_allele, '/', ALT_allele, sep = '')
            
            # Define the ith ODDS RATIO
            ith_variant_OR <- snp_feature_df[ith_relationship, 'OR.R.']
            
            # Define the ith SV FDR
            ith_FDR <- snp_feature_df[ith_relationship, 'FDR']
            
            # Define the ith plot main title (just the SV ID)
            ith_plot_title <- paste(SV_to_plot, sep = ' ')
            
              

            
            
            # LOOP OVER ETHNIC GROUPS FOR THE ITH SV
            
            # Start ethnic group loop counter
            group_loop_counter <- 0
            
            for (group_i in ethnic_groups) {
              
              
                # Update the loop counter
                group_loop_counter <- group_loop_counter + 1
              
                # Define the genotypes to use, changing with each ethnic group
                if (group_i == 'AFR') {
                      group_i_genotypes <- AFR_binary_genotypes
                }
            
                if (group_i == 'AMR') {
                      group_i_genotypes <- AMR_binary_genotypes
                }
            
                if (group_i == 'EUR') {
                      group_i_genotypes <- EUR_binary_genotypes
                }
            
                if (group_i == 'EAS') {
                      group_i_genotypes <- EAS_binary_genotypes
                }
            
                if (group_i == 'SAS') {
                      group_i_genotypes <- SAS_binary_genotypes
                }
            
                # Check if the current ethnic group has variation in the provided SV and calculate frequencies, otherwise skip to the next ethnic group
                if ( SV_to_plot %in% colnames(group_i_genotypes) ) {
                  
                  
                      # Extract genotypes for the specified SV
                      snp_i_genotypes <- group_i_genotypes[, SV_to_plot]
                      
                      # Extract # of L1 + Alu global singletons
                      snp_i_insertions <- sample_info[rownames(group_i_genotypes), singleton_column_name]
                      
                      # Define singleton frequencies across genotypes (in genotype 0/1/2, frequency = # with singleton / total samples for that genotype)
                      insertion_frequencies_genotype_0 <- sum(snp_i_insertions[snp_i_genotypes == 0] > 0) / length(snp_i_insertions[snp_i_genotypes == 0])
                      insertion_frequencies_genotype_1 <- sum(snp_i_insertions[snp_i_genotypes == 1] > 0) / length(snp_i_insertions[snp_i_genotypes == 1])
                      insertion_frequencies_genotype_2 <- sum(snp_i_insertions[snp_i_genotypes == 2] > 0) / length(snp_i_insertions[snp_i_genotypes == 2])
                      
                      # Combine frequencies into one vector
                      all.insertion.frequencies <- c(insertion_frequencies_genotype_0, insertion_frequencies_genotype_1, insertion_frequencies_genotype_2)
                      
                      # Combine genotypes and frequencies into one dataframe
                      genotypes_phenotypes <- data.frame(genotypes = c(allele_combo_1, allele_combo_2, allele_combo_3),
                                                         Insertion_Frequencies = all.insertion.frequencies,
                                                         Ancestry = group_i)
                      
                      # Start or extend df with insertion frequencies across ancestries, depending on whether this is the first group
                      if (group_loop_counter == 1) {
                          ALL.Groups.genotypes.phenotypes <- genotypes_phenotypes
                      } else {
                          ALL.Groups.genotypes.phenotypes <- rbind(ALL.Groups.genotypes.phenotypes, genotypes_phenotypes)
                      }
                        
                } else { # If the ith variants isn't in group_i
                  
                      # Fill in Insertion Frequency df with
                      genotypes_phenotypes <- data.frame(genotypes = c(allele_combo_1, allele_combo_2, allele_combo_3),
                                                         Insertion_Frequencies = NA,
                                                         Ancestry = group_i)
                      
                      # Start or extend df with insertion frequencies across ancestries, depending on whether this is the first group
                      if (group_loop_counter == 1) {
                          ALL.Groups.genotypes.phenotypes <- genotypes_phenotypes
                      } else {
                          ALL.Groups.genotypes.phenotypes <- rbind(ALL.Groups.genotypes.phenotypes, genotypes_phenotypes)
                      }
                    
                      # Skip to the next ETHNIC GROUP ITERATION
                      next
                      
                } # CLOSE IF LOOP CHECKING IF SV IS PRESENT IN THE POPULATION
              
            } # CLOSE FOR LOOP OVER ETHNIC GROUPS
            
            
            
            
            
            # CONTINUE WITH THE ITH SV PLOT
            
            # Add colors to the final results table
            ALL.Groups.genotypes.phenotypes$colors <- ancestral.group.colors
            
            # Change NaN to NA (since beeswarm can't interpret NaN)
            ALL.Groups.genotypes.phenotypes[which(ALL.Groups.genotypes.phenotypes$Insertion_Frequencies == 'NaN'), 'Insertion_Frequencies'] <- NA
            
            # Convert genotype to factor in order to preserve order in the plot
            ALL.Groups.genotypes.phenotypes$genotypes <- factor(ALL.Groups.genotypes.phenotypes$genotypes, levels = c(allele_combo_1, allele_combo_2, allele_combo_3))
          
            # Beeswarm plot
            beeswarm(ALL.Groups.genotypes.phenotypes$Insertion_Frequencies ~ ALL.Groups.genotypes.phenotypes$genotypes,
                     pch = 19, 
                     pwcol = ALL.Groups.genotypes.phenotypes$colors,
                     main = ith_plot_title,
                     xlab = c('Genotype'),
                     ylab = c('Fraction with L1/Alu Global Singletons'),
                     ylim = c(0, 1),
                     cex = 1.5
                     )
            
            # Specify axis labels
            axis(1, at = c(1, 2, 3), labels = c(allele_combo_1, allele_combo_2, allele_combo_3))
            
            # Add FDR and Odds Ratio to plot
            text(1, 1.0, paste('FDR = ', signif(ith_FDR, 3), sep = ''), cex = 1)
            text(0.79, 0.93, paste('OR = ', signif(ith_variant_OR, 3), sep = ''), cex = 1)
            
            # Legend
            legend(legend_position, 
                   legend = unique(ALL.Groups.genotypes.phenotypes$Ancestry),
                   col = unique(ALL.Groups.genotypes.phenotypes$colors), 
                   pch = 19,
                   cex = 1,
                   bty = 'n')
      
          
        } # END FOR LOOP OVER THE RELATIONSHIPS YOU WOULD LIKE TO PLOT
      
        
    # Close the plot object
    dev.off()
    
    # Nothing needs to be returned
    return(NULL)
        
  
} # END FUNCTION

plot_genotype_vs_Singleton_Freq <- function(snp_feature_df, sample_info, singleton_column_name, snp.mapping.info, dir.output, output_name, AFR_binary_genotypes, AMR_binary_genotypes, EUR_binary_genotypes, EAS_binary_genotypes, SAS_binary_genotypes, legend_position){
  
    # GENERATE 3X3 PLOTS FOR 9 SNP VS INSERTION FREQUENCY RELATIONSHIPS
  
    # SNP_FEATURE_DF: MAKE SURE ITS ONLY 5 COLUMNS, 1 IS THE SNP, 2 IS THE FEATURE (GENE,TE, ETC.), 3 IS A DESCRIPTOR OF THEIR RELATIONSHIP (OVERLAP, NEAR, ETC), 4 is the FDR, and 5 is the Odds Ratio
    # SAMPLE_INFO: Contains the number of global singletons in the column specificied by singleton_column_name
    # binary_genotype objects are from BEDMATRIX

  
  
  
    # DEFINE GENERAL PARAMETERS
  
    # Update the colnames for snp_feature_df
    colnames(snp_feature_df) <- c('SNP', 'Feature', 'Relationship', 'FDR', 'OR.R.')

    # Define vector of ethnic groups to loop over
    ethnic_groups <- c('AFR', 'AMR', 'EUR', 'EAS', 'SAS')
    
    # Define colors for points
    ancestral.group.colors <- c(rep('#648FFF', 3),
                                rep('#785EF0', 3),
                                rep('#FE6100', 3),
                                rep('#DC267F', 3),
                                rep('#FFB000', 3)
                                )
    
    
  
    
    # START THE PLOT OBJECT
    pdf(paste(dir.output, output_name, ".pdf", sep=""), width = 10, height = 10)
    par(mfrow = c(3,3), pty ="m") # 3 rows by 3 columns
    
        # Loop over the entries you would like to plot
        for (ith_relationship in 1:nrow(snp_feature_df)) {
          
            # DEFINE VARIABLES FOR THE ITH SNP
          
            # Define the ith snp to plot
            snp_to_plot <- snp_feature_df[ith_relationship, 'SNP']
            
                # Define the allele combinations for that snp
                REF_allele <- snp.mapping.info[snp_to_plot, 'REF']
                ALT_allele <- snp.mapping.info[snp_to_plot, 'ALT']
                allele_combo_1 <- paste(REF_allele, '/', REF_allele, sep = '')
                allele_combo_2 <- paste(REF_allele, '/', ALT_allele, sep = '')
                allele_combo_3 <- paste(ALT_allele, '/', ALT_allele, sep = '')
            
            # Define the ith feature
            snp_feature <- snp_feature_df[ith_relationship, 'Feature']
            
            # Define the ith snp-gene relationship
            snp_gene_relationship <- snp_feature_df[ith_relationship, 'Relationship']
            
            # Define the ith snp FDR
            ith_FDR <- snp_feature_df[ith_relationship, 'FDR']
            
            # Define the ith ODDS RATIO
            ith_variant_OR <- snp_feature_df[ith_relationship, 'OR.R.']
            
            # Define the ith plot main title
            ith_plot_title <- paste(snp_to_plot, snp_gene_relationship, snp_feature, sep = ' ')
            
            # Italicize gene 
            ith_plot_title <- substitute(paste(bold(snp_to_plot)~bold(snp_gene_relationship)~bolditalic(snp_feature)))
              

            
            
            
            # LOOP OVER ETHNIC GROUPS FOR THE ITH SNP
            
            # Start ethnic group loop counter
            group_loop_counter <- 0
            
            for (group_i in ethnic_groups) {
              
              
                # Update the loop counter
                group_loop_counter <- group_loop_counter + 1
              
                # Define the genotypes to use, changing with each ethnic group
                if (group_i == 'AFR') {
                      group_i_genotypes <- AFR_binary_genotypes
                }
            
                if (group_i == 'AMR') {
                      group_i_genotypes <- AMR_binary_genotypes
                }
            
                if (group_i == 'EUR') {
                      group_i_genotypes <- EUR_binary_genotypes
                }
            
                if (group_i == 'EAS') {
                      group_i_genotypes <- EAS_binary_genotypes
                }
            
                if (group_i == 'SAS') {
                      group_i_genotypes <- SAS_binary_genotypes
                }
            
                # Check if the current ethnic group has variation in the provided snp and calculate frequencies, otherwise skip to the next ethnic group
                if ( snp_to_plot %in% colnames(group_i_genotypes) ) {
                  
                  
                      # Extract genotypes for the specified snp
                      snp_i_genotypes <- group_i_genotypes[, snp_to_plot]
                      
                      # Extract # of L1 + Alu global singletons
                      snp_i_insertions <- sample_info[rownames(group_i_genotypes), singleton_column_name]
                      
                      # Define singleton frequencies across genotypes (in genotype 0/1/2, frequency = # with singleton / total samples)
                      insertion_frequencies_genotype_0 <- sum(snp_i_insertions[snp_i_genotypes == 0] > 0) / length(snp_i_insertions[snp_i_genotypes == 0])
                      insertion_frequencies_genotype_1 <- sum(snp_i_insertions[snp_i_genotypes == 1] > 0) / length(snp_i_insertions[snp_i_genotypes == 1])
                      insertion_frequencies_genotype_2 <- sum(snp_i_insertions[snp_i_genotypes == 2] > 0) / length(snp_i_insertions[snp_i_genotypes == 2])
                      
                      # Combine frequencies into one vector
                      all.insertion.frequencies <- c(insertion_frequencies_genotype_0, insertion_frequencies_genotype_1, insertion_frequencies_genotype_2)
                      
                      # Combine genotypes and frequencies into one dataframe
                      genotypes_phenotypes <- data.frame(genotypes = c(allele_combo_1, allele_combo_2, allele_combo_3),
                                                         Insertion_Frequencies = all.insertion.frequencies,
                                                         Ancestry = group_i)
                      
                      # Start or extend df with insertion frequencies across ancestries, depending on whether this is the first group
                      if (group_loop_counter == 1) {
                          ALL.Groups.genotypes.phenotypes <- genotypes_phenotypes
                      } else {
                          ALL.Groups.genotypes.phenotypes <- rbind(ALL.Groups.genotypes.phenotypes, genotypes_phenotypes)
                      }
                        
                } else { # If the ith snp isn't in group_i
                  
                      # Fill in Insertion Frequency df with
                      genotypes_phenotypes <- data.frame(genotypes = c(allele_combo_1, allele_combo_2, allele_combo_3),
                                                         Insertion_Frequencies = NA,
                                                         Ancestry = group_i)
                      
                      # Start or extend df with insertion frequencies across ancestries, depending on whether this is the first group
                      if (group_loop_counter == 1) {
                          ALL.Groups.genotypes.phenotypes <- genotypes_phenotypes
                      } else {
                          ALL.Groups.genotypes.phenotypes <- rbind(ALL.Groups.genotypes.phenotypes, genotypes_phenotypes)
                      }
                    
                      # Skip to the next ETHNIC GROUP ITERATION
                      next
                      
                } # CLOSE IF LOOP CHECKING IF SNP IS PRESENT IN THE POPULATION
              
            } # CLOSE FOR LOOP OVER ETHNIC GROUPS
            
            
            
            
            
            # CONTINUE WITH THE ITH SNP PLOT
            
            # Add colors to the final results table
            ALL.Groups.genotypes.phenotypes$colors <- ancestral.group.colors
            
            # Change NaN to NA (since beeswarm can't interpret NaN)
            ALL.Groups.genotypes.phenotypes[which(ALL.Groups.genotypes.phenotypes$Insertion_Frequencies == 'NaN'), 'Insertion_Frequencies'] <- NA
            
            # Convert genotype to factor in order to preserve order in the plot
            ALL.Groups.genotypes.phenotypes$genotypes <- factor(ALL.Groups.genotypes.phenotypes$genotypes, levels = c(allele_combo_1, allele_combo_2, allele_combo_3))
          
            # Beeswarm plot
            beeswarm(ALL.Groups.genotypes.phenotypes$Insertion_Frequencies ~ ALL.Groups.genotypes.phenotypes$genotypes,
                     pch = 19, 
                     pwcol = ALL.Groups.genotypes.phenotypes$colors,
                     main = ith_plot_title,
                     xlab = c('Genotype'),
                     ylab = c('Fraction with L1/Alu Global Singletons'),
                     ylim = c(0, 1),
                     cex = 1.5
                     )
            
            # Specify axis labels
            axis(1, at = c(1, 2, 3), labels = c(allele_combo_1, allele_combo_2, allele_combo_3))
            
            # Add FDR and Odds Ratio to plot
            text(1, 1.0, paste('FDR = ', signif(ith_FDR, 3), sep = ''), cex = 1)
            text(0.79, 0.93, paste('OR = ', signif(ith_variant_OR, 3), sep = ''), cex = 1)
            
            # Legend
            legend(legend_position, 
                   legend = unique(ALL.Groups.genotypes.phenotypes$Ancestry),
                   col = unique(ALL.Groups.genotypes.phenotypes$colors), 
                   pch = 19,
                   cex = 1,
                   bty = 'n')
      
          
        } # END FOR LOOP OVER THE RELATIONSHIPS YOU WOULD LIKE TO PLOT
      
        
    # Close the plot object
    dev.off()
    
    # Nothing needs to be returned
    return(NULL)
        
  
} # END FUNCTION
