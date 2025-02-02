##########################################################################
# This R script summarizes the development of all 04 - 13 prefixed scripts
# Particularly it only uses 09 and 13.
#
# Created by: Roman Cerny
# Date:       2021-07-12
#
##########################################################################

library(dplyr)
library(stringdist)
library(plotly)

require(parallel)
require(pbapply)
require(emoa)
require(mco)
require(ggplot2)

library(BiocManager)
library(Biostrings)

# Specifying files to work with
#
# e.g. list(
#         c(<1_organism_name>, <1_organism_predictions_path>, <1_peptides_excl_organism_path>),
#         c(<2_organism_name>, <2_organism_predictions_path>, <2_peptides_excl_organism_path>),
#         c(<3_organism_name>, <3_organism_predictions_path>, <3_peptides_excl_organism_path>)
#      )
#
###############################

organisms_data_paths <- list(
  c('01_EBV', './data/data sets/01_EBV/pred_peptides.rds',
    './data/data sets/01_EBV/epitopes_except_organism.rds'),
    c('02_HepC', './data/data sets/02_HepC/pred_peptides.rds',
      './data/data sets/02_HepC/epitopes_except_organism.rds'),
    c('03_Ovolvulus', './data/data sets/03_Ovolvulus/pred_peptides.rds',
      './data/data sets/03_Ovolvulus/epitopes_except_organism.rds')
)

peptides_path <- './data/data sets/peptides.rds'

# Algorithm params
##################

substitution_matrix <- 'PAM30'
pairwise_align_type <- 'local' # Smith-Waterman
gap_opening <- 5
gap_extension <- 2

# Calculate all peptides self-alignment score if not present
############################################################

peptides <- readRDS(file=peptides_path)

if (!"Self_alignment_score" %in% colnames(peptides)) {
  peptides$Self_alignment_score <- pairwiseAlignment(unlist(peptides$Info_peptide),
                                                     unlist(peptides$Info_peptide),
                                                     type=pairwise_align_type,
                                                     substitutionMatrix=substitution_matrix, 
                                                     gapOpening=gap_opening, 
                                                     gapExtension=gap_extension,
                                                     scoreOnly = TRUE)
  
  saveRDS(peptides, file=peptides_path)
}


for (organism_data_path in organisms_data_paths) {
  
  organism_predictions <- readRDS(file=organism_data_path[2])
  peptides_excl_organism <- readRDS(file=organism_data_path[3])
  
  # Calculate organism self-alignment score if not present
  ############################################################
  
  if (!"Self_alignment_score" %in% colnames(organism_predictions)) {
    organism_predictions$Self_alignment_score <- pairwiseAlignment(unlist(organism_predictions$Sequence),
                                                                   unlist(organism_predictions$Sequence),
                                                                   type=pairwise_align_type,
                                                                   substitutionMatrix=substitution_matrix, 
                                                                   gapOpening=gap_opening, 
                                                                   gapExtension=gap_extension,
                                                                   scoreOnly = TRUE)
    
    saveRDS(organism_predictions, file=organism_data_path[2])
  }
  
  
  # Normalize column names
  colnames(organism_predictions) <- c('Protein', 
                                      'Start_pos', 
                                      'End_pos', 
                                      'Length', 
                                      'Probability', 
                                      'Sequence', 
                                      'CandidateID',
                                      'Self_alignment_score')
  
  # Calculating alignment scores
  ##############################
  
  cores <- detectCores()
  
  cl <- makeCluster(cores-1)
  parallel::clusterExport(cl=cl, varlist=c('organism_predictions', 
                                           'peptides_excl_organism',
                                           'pairwise_align_type',
                                           'gap_opening',
                                           'gap_extension',
                                           'substitution_matrix',
                                           'peptides'))
  
  parallel::clusterEvalQ(cl=cl, library(Biostrings))
  
  system.time({
    
    result = pblapply(cl=cl,
                      X=1:nrow(organism_predictions),
                      FUN=function(idx) {
                        protein <- organism_predictions[idx, 'Protein']
                        start_pos <- organism_predictions[idx, 'Start_pos']
                        end_pos <- organism_predictions[idx, 'End_pos']
                        pep_length <- organism_predictions[idx, 'Length']
                        probability <- organism_predictions[idx, 'Probability']
                        sequence <- organism_predictions[idx, 'Sequence']
                        self_alignment_score_pred <- organism_predictions[idx, 'Self_alignment_score']
                        
                        result_all <- as.data.frame(matrix(nrow=nrow(peptides_excl_organism), ncol=9))
                        
                        for (peptide_idx in 1:nrow(peptides_excl_organism)) {
                          info_PepID <- peptides_excl_organism[peptide_idx, 'Info_PepID']
                          info_peptide <- peptides_excl_organism[peptide_idx, 'Info_peptide']
                          self_alignment_score_peptide <- peptides[which(peptides$Info_PepID==unlist(info_PepID)), 'Self_alignment_score']
                          
                          alignment_score <- pairwiseAlignment(unlist(sequence), unlist(info_peptide),
                                                     type=pairwise_align_type,
                                                     substitutionMatrix=substitution_matrix, 
                                                     gapOpening=gap_opening, 
                                                     gapExtension=gap_extension,
                                                     scoreOnly = TRUE)
                          
                          uniqueness_score <- 1 - (alignment_score / min(unlist(self_alignment_score_peptide),
                                                                         unlist(self_alignment_score_pred)))
                          
                          result_all[peptide_idx, 1] <- protein
                          result_all[peptide_idx, 2] <- start_pos
                          result_all[peptide_idx, 3] <- end_pos
                          result_all[peptide_idx, 4] <- pep_length
                          result_all[peptide_idx, 5] <- probability
                          result_all[peptide_idx, 6] <- sequence
                          result_all[peptide_idx, 7] <- uniqueness_score
                          result_all[peptide_idx, 8] <- alignment_score
                          result_all[peptide_idx, 9] <- info_peptide
                          result_all[peptide_idx, 10] <- info_PepID
                        }
                        
                        min_uniqueness <- 999999999
                        min_uniqueness_row_idx <- 0
                        
                        # Find min uniqueness score
                        for (peptide_idx in 1:nrow(result_all)) {
                          if (result_all[peptide_idx, 7] < min_uniqueness) {
                            min_uniqueness <- result_all[peptide_idx, 7]
                            min_uniqueness_row_idx <- peptide_idx
                          }
                        }
                        
                        # get only 1 shortest distance for each protein
                        result = as.data.frame(matrix(nrow=1, ncol=9))
                        result[1, 1] <- result_all[min_uniqueness_row_idx, 1]
                        result[1, 2] <- result_all[min_uniqueness_row_idx, 2]
                        result[1, 3] <- result_all[min_uniqueness_row_idx, 3]
                        result[1, 4] <- result_all[min_uniqueness_row_idx, 4]
                        result[1, 5] <- result_all[min_uniqueness_row_idx, 5]
                        result[1, 6] <- result_all[min_uniqueness_row_idx, 6]
                        result[1, 7] <- result_all[min_uniqueness_row_idx, 7]
                        result[1, 8] <- result_all[min_uniqueness_row_idx, 8]
                        result[1, 9] <- result_all[min_uniqueness_row_idx, 9]
                        result[1, 10] <- result_all[min_uniqueness_row_idx, 10]
                        
                        return(result)
                      })
    
    result <- dplyr::bind_rows(result)
    stopCluster(cl)
    
    colnames(result) <- c('Protein',
                          'Start_pos',
                          'End_pos',
                          'Length',
                          'Probability',
                          'Sequence',
                          'Uniqueness_score',
                          'Alignment_score',
                          'Info_peptide',
                          'Info_PepID')
    
    print('Writing result to the file...')
    
    output_file_name <- paste('./output/predictions_pairwiseAlignment_', 
                              organism_data_path[1], '_', 
                              substitution_matrix, '_', 
                              pairwise_align_type, '_', 
                              Sys.time(), sep='')
    output_file_name <- gsub(' ', '_', output_file_name, fixed=TRUE)
    output_file_name <- gsub(':', '_', output_file_name, fixed=TRUE)
    
    write.csv(result,
              paste(output_file_name, '.csv', sep=''), 
              row.names=FALSE)
    
    ## persist histogram in file
    #png(paste(output_file_name, '.png', sep=''))
    #hist(result$Uniqueness_score)
    #dev.off()
  })
  
  
  # Non-dominated ranking & plotting
  ##################################
  
  PLOT_SEQUENCE <- FALSE
  file_name <- paste(output_file_name, '.csv', sep='')
  
  df <- read.csv(file=file_name)
  
  df_complemented <- df %>%
    mutate(probability_complement = 1 - Probability) %>%
    mutate(uniqueness_score_complement = 1 - Uniqueness_score) %>%
    select(Protein, Sequence, 
           Probability, Uniqueness_score, 
           probability_complement, 
           uniqueness_score_complement, 
           Info_peptide, Info_PepID)
  
  df_objectives <- df_complemented %>%
    select(probability_complement, uniqueness_score_complement)
  
  m_objectives <- t(as.matrix(df_objectives))
  
  nondominated <- nondominated_points(m_objectives)
  
  rank <- nds_rank(m_objectives)
  
  file_name_split <- unlist(strsplit(file_name, '_'))
  
  if (PLOT_SEQUENCE) {
    pl <- ggplot(data=df_complemented, aes(x=Probability, y=Uniqueness_score)) +
      ggtitle(paste(organism_data_path[1], substitution_matrix, sep=' ')) +
      geom_point(aes(text=sprintf("Protein ID: %s Sequence: %s<br>Info_PepID: %s Info_peptide: %s", 
                                  Protein, Sequence, 
                                  Info_PepID, Info_peptide)))
  } else {
    pl <- ggplot(data=df_complemented, aes(x=Probability, y=Uniqueness_score)) +
      ggtitle(paste(organism_data_path[1], substitution_matrix, sep=' ')) +
      geom_point(aes(text=sprintf("Protein ID: %s <br>Info_PepID: %s", 
                                  Protein, 
                                  Info_PepID)))
  }
  
  # Plot Pareto fronts
  for (idx in 1:max(rank)) {
    color <- 'blue'
    if (!idx %% 2) {
      color <- 'green'
    }
    
    # Prepare Pareto points
    pareto_front_points_indicies <- which(rank %in% c(idx))
    df_pareto_front_points <- df_objectives[pareto_front_points_indicies, 
                                            c('probability_complement', 'uniqueness_score_complement')]
    
    
    df_pareto_front_points_not_inverted <- df_pareto_front_points %>%
      mutate(Probability = 1 - probability_complement) %>%
      mutate(Uniqueness_score = 1 - uniqueness_score_complement) %>%
      select(Probability, Uniqueness_score)
    
    pl <- pl + geom_line(data=df_pareto_front_points_not_inverted, 
                         aes(x=Probability, y=Uniqueness_score), 
                         color=color)
  }
  
  chart <- ggplotly(pl) 
  
  # Save plot in interactive file
  htmlwidgets::saveWidget(chart, selfcontained=TRUE, paste(output_file_name, '.html', sep=''))
}
