library(dplyr)
library(stringdist)
require(parallel)
require(pbapply)

Ov_predictions <- read.csv(file = './data/Ov_predictions.csv')
head(Ov_predictions)

peptides_no_Ov <- read.csv(file = './data/peptides_no_Ov.csv')
head(peptides_no_Ov)

lv_del_penalty_w <- 0.5
lv_ins_penalty_w <- 1
lv_sub_penalty_w <- 1

# https://www.r-bloggers.com/2017/05/5-ways-to-measure-running-time-of-r-code/
  
cores <- detectCores()
cl <- makeCluster(cores-1)
parallel::clusterExport(cl= cl, varlist = c("Ov_predictions", 
                                            "peptides_no_Ov",
                                            "lv_del_penalty_w",
                                            "lv_ins_penalty_w",
                                            "lv_sub_penalty_w"))

parallel::clusterEvalQ(cl= cl, library(stringdist))

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
                      
                      result = as.data.frame(matrix(nrow=nrow(peptides_no_Ov), ncol=10))
                      
                      for (peptide_idx in 1:nrow(peptides_no_Ov)) {
                        info_PepID <- peptides_no_Ov[peptide_idx, 'Info_PepID']
                        info_peptide <- peptides_no_Ov[peptide_idx, 'Info_peptide']
                        
                        #lev_distances <- data.frame()
                        lev_distances <- list()
                        
                        a_to_b = stringdist(sequence,
                                            info_peptide,
                                            method = c("lv"),
                                            useBytes = FALSE,
                                            weight = c(d = lv_del_penalty_w, i = lv_ins_penalty_w, s = lv_sub_penalty_w),
                                            nthread = getOption("sd_num_thread"))
                        
                        b_to_a = stringdist(info_peptide,
                                            sequence,
                                            method = c("lv"),
                                            useBytes = FALSE,
                                            weight = c(d = lv_del_penalty_w, i = lv_ins_penalty_w, s = lv_sub_penalty_w),
                                            nthread = getOption("sd_num_thread"))
                        
                        result[peptide_idx, 1] <- protein
                        result[peptide_idx, 2] <- start_pos
                        result[peptide_idx, 3] <- end_pos
                        result[peptide_idx, 4] <- pep_length
                        result[peptide_idx, 5] <- probability
                        result[peptide_idx, 6] <- sequence
                        result[peptide_idx, 7] <- a_to_b
                        result[peptide_idx, 8] <- b_to_a
                        result[peptide_idx, 9] <- info_peptide
                        result[peptide_idx, 10] <- info_PepID
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
                        "A_to_B_lv_dist",
                        "B_to_A_lv_dist",
                        "Info_peptide",
                        "Info_PepID")

  print("Writing result to the file...")
    
  write.csv(result,
            "./output/predictions_lev_dist.csv", 
            row.names = FALSE)
})

View(result)

res <- read.csv(file = './output/predictions_lev_dist.csv')

hist(res$A_to_B_lv_dist)

View(res)







