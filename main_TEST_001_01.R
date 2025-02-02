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
    './data/data sets/01_EBV/epitopes_except_organism.rds')
  #  c('02_HepC', './data/data sets/02_HepC/pred_peptides.rds',
  #    './data/data sets/02_HepC/epitopes_except_organism.rds'),
  #  c('03_Ovolvulus', './data/data sets/03_Ovolvulus/pred_peptides.rds',
  #    './data/data sets/03_Ovolvulus/epitopes_except_organism.rds')
)

peptides_path <- './data/data sets/peptides.rds'

# Algorithm params
##################

substitution_matrix <- 'PAM30'
pairwise_align_type <- 'local' # 'overlap'
gap_opening <- 5
gap_extension <- 2

# Calculate all peptides self-alignment score if not present
############################################################

peptides <- readRDS(file=peptides_path)

if (!"Self_alignment_score" %in% colnames(peptides)) {
  peptides$Self_alignment_score <- pairwiseAlignment(peptides$Info_peptide,
                                                     peptides$Info_peptide,
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
  
  #organism_predictions <- organism_predictions[1:50,] # 53
  peptides_excl_organism <- peptides_excl_organism[1:10,] # 13000
  
  
  
  
  # Calculate organism self-alignment score if not present
  ############################################################
  
  if (!"Self_alignment_score" %in% colnames(organism_predictions)) {
    organism_predictions$Self_alignment_score <- pairwiseAlignment(organism_predictions$Sequence,
                                                                   organism_predictions$Sequence,
                                                                   type=pairwise_align_type,
                                                                   substitutionMatrix=substitution_matrix, 
                                                                   gapOpening=gap_opening, 
                                                                   gapExtension=gap_extension,
                                                                   scoreOnly = TRUE)
    
    ### !!!!!!!!!!!!!!!!!!!!!!!! %>% head()%>% head()%>% head() ###  ##saveRDS(organism_predictions, file=organism_data_path[2])
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
  
  #cores <- detectCores()
  
  #cl <- makeCluster(cores-1)
  #parallel::clusterExport(cl=cl, varlist=c('organism_predictions', 
  #                                         'peptides_excl_organism',
  #                                         'pairwise_align_type',
  #                                         'gap_opening',
  #                                         'gap_extension',
  #                                         'substitution_matrix',
  #                                         'peptides'))
  
  #parallel::clusterEvalQ(cl=cl, library(Biostrings))
  
  #system.time({
    
  #  result = pblapply(cl=cl,
  #                    X=1:nrow(organism_predictions),
  #                    FUN=function(idx) {
  idx<-1
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

                          print('---------------')
                          print(self_alignment_score_peptide)
                          print('---------------')
                          
                          
                          score <- pairwiseAlignment('GKAQPIQSSHLSSMSPTQPISHEEQPRYEDPDAPLDLSLHPDVAAQPAPQAPYQGYQEPPAPQAPYQGYQEPPAPQAPYQGYQEPPAPQAPYQGYQEPPAPQAPYQGYQEPPAPQAPYQGYQEPPAPQAPYQGYQEPPAPGLQSSSYPGYAGPWTPRSQHPCYRHPWAPWSQDPVHGHTQGPWDPRAPHLPPQWDGSAGHGQDQVSQFPHLQSETGPPRLQLSLVPLVSSSAPSWSSPQPRAPIRPIPT', 'APQPGPQPPQPPQPQPEAPAPQPPAG',
                                                     type=pairwise_align_type,
                                                     substitutionMatrix=substitution_matrix, 
                                                     gapOpening=gap_opening, 
                                                     gapExtension=gap_extension,
                                                     scoreOnly = TRUE)

                          uniqueness <- 1 - (score / min(
                                                     unlist(self_alignment_score_peptide),
                                                     unlist(self_alignment_score_pred)))
                          
                          print(score)
                          print(min(peptides$Self_alignment_score))
                          print(min(unlist(self_alignment_score)))
                          print(min(min(peptides$Self_alignment_score), unlist(self_alignment_score)))
                          print(uniqueness)
                          print('*************************************')
                          
                          
                          
                                                    result_all[peptide_idx, 1] <- protein
                          result_all[peptide_idx, 2] <- start_pos
                          result_all[peptide_idx, 3] <- end_pos
                          result_all[peptide_idx, 4] <- pep_length
                          result_all[peptide_idx, 5] <- probability
                          result_all[peptide_idx, 6] <- sequence
                          result_all[peptide_idx, 7] <- uniqueness
                          result_all[peptide_idx, 8] <- info_peptide
                          result_all[peptide_idx, 9] <- info_PepID
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
                        
#                        return(result)
#                      })
    
#    result <- dplyr::bind_rows(result)
#    stopCluster(cl)
    
##    colnames(result) <- c('Protein',
##                          'Start_pos',
##                          'End_pos',
##                          'Length',
##                          'Probability',
##                          'Sequence',
##                          'Uniqueness_score',
##                          'Info_peptide',
##                          'Info_PepID')
##    
##    print('Writing result to the file...')
##    
##    output_file_name <- paste('./output/predictions_pairwiseAlignment_', 
##                              organism_data_path[1], '_', 
##                              substitution_matrix, '_', 
##                              pairwise_align_type, '_', 
##                              Sys.time(), sep='')
##    output_file_name <- gsub(' ', '_', output_file_name, fixed=TRUE)
##    output_file_name <- gsub(':', '_', output_file_name, fixed=TRUE)
##    
##    write.csv(result,
##              paste(output_file_name, '.csv', sep=''), 
##              row.names=FALSE)
##    
    ## persist histogram in file
    #png(paste(output_file_name, '.png', sep=''))
    #hist(result$Uniqueness_score)
    #dev.off()
#  })

}

