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

substitution_matrix <- 'PAM30'
#substitution_matrix <- 'BLOSUM45'

gap_opening <- 3
gap_extension <- 1

cores <- detectCores()
cl <- makeCluster(cores-1)
parallel::clusterExport(cl= cl, varlist = c("Ov_predictions", 
                                            "peptides_no_Ov",
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
                      
                      result = as.data.frame(matrix(nrow=nrow(peptides_no_Ov), ncol=9))
                      
                      for (peptide_idx in 1:nrow(peptides_no_Ov)) {
                        info_PepID <- peptides_no_Ov[peptide_idx, 'Info_PepID']
                        info_peptide <- peptides_no_Ov[peptide_idx, 'Info_peptide']
                        
                        pair_align = pairwiseAlignment(sequence, info_peptide,
                                                       type = "overlap",
                                                       substitutionMatrix = substitution_matrix, 
                                                       gapOpening = gap_opening, 
                                                       gapExtension = gap_extension)
                        
                        
                        result[peptide_idx, 1] <- protein
                        result[peptide_idx, 2] <- start_pos
                        result[peptide_idx, 3] <- end_pos
                        result[peptide_idx, 4] <- pep_length
                        result[peptide_idx, 5] <- probability
                        result[peptide_idx, 6] <- sequence
                        result[peptide_idx, 7] <- pair_align@score
                        result[peptide_idx, 8] <- info_peptide
                        result[peptide_idx, 9] <- info_PepID
                      }
                      
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
  
  output_file_name <- paste('./output/predictions_smith_waterman_', substitution_matrix, '_', Sys.time(), '.csv', sep="")
  output_file_name <- gsub(' ', '_', output_file_name, fixed = TRUE)
  output_file_name <- gsub(':', '_', output_file_name, fixed = TRUE)
  
  write.csv(result,
            output_file_name, 
            row.names = FALSE)
})

View(result)

#res <- read.csv(file = output_file_name)

hist(result$pair_align_score)

#View(res)

