library(dplyr)
library(stringdist)
require(parallel)
require(pbapply)

library(BiocManager)
library(Biostrings)

Ov_predictions <- read.csv(file = './data/Ov_predictions.csv')
head(Ov_predictions)

peptides_no_Ov <- read.csv(file = './data/peptides_no_Ov.csv')
head(peptides_no_Ov)


# ,'BLOSUM45'
# ,'BLOSUM50'
# ,'BLOSUM62'
# ,'BLOSUM80'
# ,'BLOSUM100'
# 
# ,'PAM30'
# ,'PAM40'
# ,'PAM70'
# ,'PAM120' 
# ,'PAM250' 

#substitution_matrices <- c('PAM30', 'BLOSUM45')
substitution_matrices <- c('PAM30')
#substitution_matrices <- c('PAM30','BLOSUM45' ,'PAM40','BLOSUM50' ,'PAM70','BLOSUM62' ,'PAM120','BLOSUM80' ,'PAM250','BLOSUM100')

pairwise_align_type <- "overlap"
#pairwise_align_type <- "local"
gap_opening <- 5
gap_extension <- 2

cores <- detectCores()

for (substitution_matrix in substitution_matrices) {
  cl <- makeCluster(cores-1)
  parallel::clusterExport(cl= cl, varlist = c("Ov_predictions", 
                                              "peptides_no_Ov",
                                              "pairwise_align_type",
                                              "gap_opening",
                                              "gap_extension",
                                              "substitution_matrix"))
  
  parallel::clusterEvalQ(cl= cl, library(Biostrings))
  
  system.time({
    
    result = pblapply(cl = cl,
                      X = 1:nrow(Ov_predictions),
                      FUN = function(idx) {
                        protein <- Ov_predictions[idx, 'Protein']
                        start_pos <- Ov_predictions[idx, 'Start_pos']
                        end_pos <- Ov_predictions[idx, 'End_pos']
                        pep_length <- Ov_predictions[idx, 'Length']
                        probability <- Ov_predictions[idx, 'Probability']
                        sequence <- Ov_predictions[idx, 'Sequence']
                        
                        result_all = as.data.frame(matrix(nrow=nrow(peptides_no_Ov), ncol=9))
                        
                        max_similarity_score <- 0
                        len_of_max_similarity_score_peptide_string <- 0
                        
                        for (peptide_idx in 1:nrow(peptides_no_Ov)) {
                          info_PepID <- peptides_no_Ov[peptide_idx, 'Info_PepID']
                          info_peptide <- peptides_no_Ov[peptide_idx, 'Info_peptide']
                          
                          pair_align = pairwiseAlignment(sequence, info_peptide,
                                                         type = pairwise_align_type,
                                                         substitutionMatrix = substitution_matrix, 
                                                         gapOpening = gap_opening, 
                                                         gapExtension = gap_extension)
                          
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
    
    result = dplyr::bind_rows(result)
    stopCluster(cl)
    
    colnames(result) <- c("Protein",
                          "Start_pos",
                          "End_pos",
                          "Length",
                          "Probability",
                          "Sequence",
                          "pair_align_score",
                          "Info_peptide",
                          "Info_PepID")
    
    print("Writing result to the file...")
    
    output_file_name <- paste('./output/predictions_pairwiseAlignment_', substitution_matrix, '_', pairwise_align_type, '_', Sys.time(), sep='')
    output_file_name <- gsub(' ', '_', output_file_name, fixed = TRUE)
    output_file_name <- gsub(':', '_', output_file_name, fixed = TRUE)
    
    write.csv(result,
              paste(output_file_name, '.csv', sep=''), 
              row.names = FALSE)
    
    # persist histogram in file
    png(paste(output_file_name, '.png', sep=''))
    hist(result$pair_align_score)
    dev.off()
  })
}


View(result)

#res <- read.csv(file = output_file_name)

#hist(result$pair_align_score)

#View(res)

