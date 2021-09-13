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

organisms_data <- list(
  c('01_EBV', './data/data sets/01_EBV/pred_peptides.rds',
    './data/data sets/01_EBV/epitopes_except_organism.rds'),
  c('02_HepC', './data/data sets/02_HepC/pred_peptides.rds',
    './data/data sets/02_HepC/epitopes_except_organism.rds'),
  c('03_Ovolvulus', './data/data sets/03_Ovolvulus/pred_peptides.rds',
    './data/data sets/03_Ovolvulus/epitopes_except_organism.rds')
)


# Algorithm params
##################

substitution_matrices <- c('PAM30')
pairwise_align_type <- 'overlap'
gap_opening <- 5
gap_extension <- 2


for (organism_data in organisms_data) {

  # Calculating alignment scores
  ##############################
  
  organism_predictions <- readRDS(file=organism_data[2])
  peptides_excl_organism <- readRDS(file=organism_data[3])
  
  # normalize column names
  colnames(organism_predictions) <- c('Protein', 
                                      'Start_pos', 
                                      'End_pos', 
                                      'Length', 
                                      'Probability', 
                                      'Sequence', 
                                      'CandidateID')

  cores <- detectCores()
  
  for (substitution_matrix in substitution_matrices) {
    cl <- makeCluster(cores-1)
    parallel::clusterExport(cl=cl, varlist=c('organism_predictions', 
                                              'peptides_excl_organism',
                                              'pairwise_align_type',
                                              'gap_opening',
                                              'gap_extension',
                                              'substitution_matrix'))
    
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
                          
                          result_all <- as.data.frame(matrix(nrow=nrow(peptides_excl_organism), ncol=9))
                          
                          max_similarity_score <- 0
                          len_of_max_similarity_score_peptide_string <- 0
                          
                          for (peptide_idx in 1:nrow(peptides_excl_organism)) {
                            info_PepID <- peptides_excl_organism[peptide_idx, 'Info_PepID']
                            info_peptide <- peptides_excl_organism[peptide_idx, 'Info_peptide']
                            
                            pair_align = pairwiseAlignment(sequence, info_peptide,
                                                           type=pairwise_align_type,
                                                           substitutionMatrix=substitution_matrix, 
                                                           gapOpening=gap_opening, 
                                                           gapExtension=gap_extension)
                            
                            score <- pair_align@score 
                            
                            if (score > max_similarity_score) {
                              max_similarity_score <- score
                              len_of_max_similarity_score_peptide_string <- nchar(info_peptide)
                            }
                            
                            result_all[peptide_idx, 1] <- protein
                            result_all[peptide_idx, 2] <- start_pos
                            result_all[peptide_idx, 3] <- end_pos
                            result_all[peptide_idx, 4] <- pep_length
                            result_all[peptide_idx, 5] <- probability
                            result_all[peptide_idx, 6] <- sequence
                            result_all[peptide_idx, 7] <- score
                            result_all[peptide_idx, 8] <- info_peptide
                            result_all[peptide_idx, 9] <- info_PepID
                          }
                          
                          max_score <- 0
                          max_score_row_idx <- 0
                          
                          # length normalized alignment score
                          for (peptide_idx in 1:nrow(result_all)) {
                            result_all[peptide_idx, 7] <- result_all[peptide_idx, 7] / 
                              min(c(result_all[peptide_idx, 4],
                                    len_of_max_similarity_score_peptide_string))
                            
                            if (result_all[peptide_idx, 7] > max_score) {
                              max_score <- result_all[peptide_idx, 7]
                              max_score_row_idx <- peptide_idx
                            }
                          }
                          
                          # get only 1 shortest distance for each protein
                          result = as.data.frame(matrix(nrow=1, ncol=9))
                          result[1, 1] <- result_all[max_score_row_idx, 1]
                          result[1, 2] <- result_all[max_score_row_idx, 2]
                          result[1, 3] <- result_all[max_score_row_idx, 3]
                          result[1, 4] <- result_all[max_score_row_idx, 4]
                          result[1, 5] <- result_all[max_score_row_idx, 5]
                          result[1, 6] <- result_all[max_score_row_idx, 6]
                          result[1, 7] <- result_all[max_score_row_idx, 7]
                          result[1, 8] <- result_all[max_score_row_idx, 8]
                          result[1, 9] <- result_all[max_score_row_idx, 9]
                          
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
                            'pair_align_score',
                            'Info_peptide',
                            'Info_PepID')
      
      print('Writing result to the file...')
      
      output_file_name <- paste('./output/predictions_pairwiseAlignment_', 
                                organism_data[1], '_', 
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
      #hist(result$pair_align_score)
      #dev.off()
    })
  }
  
  
  # Non-dominated ranking & plotting
  ##################################
  
  PLOT_SEQUENCE <- FALSE
  file_name <- paste(output_file_name, '.csv', sep='')
  
  df <- read.csv(file=file_name)
  
  df_complemented_negated <- df %>%
    mutate(probability_complement = 1 - Probability) %>%
    select(Protein, Sequence, 
           Probability, pair_align_score, 
           probability_complement, 
           Info_peptide, Info_PepID)
  
  df_objectives <- df_complemented_negated %>%
    select(probability_complement, pair_align_score)
  
  m_objectives <- t(as.matrix(df_objectives))
  
  nondominated <- nondominated_points(m_objectives)
  
  rank <- nds_rank(m_objectives)
  
  file_name_split <- unlist(strsplit(file_name, '_'))
  
  if (PLOT_SEQUENCE) {
    pl <- ggplot(data=df_complemented_negated, aes(x=Probability, y=pair_align_score)) +
      ggtitle(paste(organism_data[1], substitution_matrix, sep=' ')) +
      geom_point(aes(text=sprintf("Protein ID: %s Sequence: %s<br>Info_PepID: %s Info_peptide: %s", 
                                  Protein, Sequence, 
                                  Info_PepID, Info_peptide)))
  } else {
    pl <- ggplot(data=df_complemented_negated, aes(x=Probability, y=pair_align_score)) +
      ggtitle(paste(organism_data[1], substitution_matrix, sep=' ')) +
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
                                            c('probability_complement', 'pair_align_score')]
    
    
    df_pareto_front_points_not_inverted <- df_pareto_front_points %>%
      mutate(Probability = 1 - probability_complement) %>%
      select(Probability, pair_align_score)
    
    pl <- pl + geom_line(data=df_pareto_front_points_not_inverted, 
                         aes(x=Probability, y=pair_align_score), 
                         color=color)
  }
  
  chart <- ggplotly(pl) 
  
  # Save plot in interactive file
  htmlwidgets::saveWidget(chart, selfcontained=TRUE, paste(output_file_name, '.html', sep=''))
}
