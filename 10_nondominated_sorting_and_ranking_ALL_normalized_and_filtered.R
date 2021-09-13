require("dplyr")
require("emoa")
require("mco")
require("ggplot2")
require("rPref")

# './data/predictions_smith_waterman_PAM30_2021-06-16_23_18_07.csv'
# './data/predictions_smith_waterman_PAM40_2021-06-17_01_40_32.csv'
# './data/predictions_smith_waterman_PAM70_2021-06-17_04_02_50.csv'
# './data/predictions_smith_waterman_PAM120_2021-06-17_06_25_28.csv'
# './data/predictions_smith_waterman_PAM250_2021-06-17_08_48_16.csv'
#
# './data/predictions_smith_waterman_BLOSUM45_2021-06-17_00_29_11.csv'
# './data/predictions_smith_waterman_BLOSUM50_2021-06-17_02_51_42.csv'
# './data/predictions_smith_waterman_BLOSUM62_2021-06-17_05_14_43.csv'
# './data/predictions_smith_waterman_BLOSUM80_2021-06-17_07_36_47.csv'
# './data/predictions_smith_waterman_BLOSUM100_2021-06-17_10_00_49.csv'

file_names <- c('./data/predictions_smith_waterman_PAM30_2021-06-16_23_18_07.csv',
                './data/predictions_smith_waterman_BLOSUM45_2021-06-17_00_29_11.csv',
                './data/predictions_smith_waterman_PAM40_2021-06-17_01_40_32.csv',
                './data/predictions_smith_waterman_BLOSUM50_2021-06-17_02_51_42.csv',
                './data/predictions_smith_waterman_PAM70_2021-06-17_04_02_50.csv',
                './data/predictions_smith_waterman_BLOSUM62_2021-06-17_05_14_43.csv',
                './data/predictions_smith_waterman_PAM120_2021-06-17_06_25_28.csv',
                './data/predictions_smith_waterman_BLOSUM80_2021-06-17_07_36_47.csv',
                './data/predictions_smith_waterman_PAM250_2021-06-17_08_48_16.csv',
                './data/predictions_smith_waterman_BLOSUM100_2021-06-17_10_00_49.csv')

png('./output/pareto_fronts_ALL.png', width=2000, height=4500, res=300)
par(mfrow=c(length(file_names)/2, 2))

for (file_name in file_names) {
  df <- read.csv(file = file_name)
  
  df_objectives <- df %>%
    mutate(probability_complement = 1 - Probability,
           distance_neg = -pair_align_score) %>%
    select(probability_complement, distance_neg) #%>% head() 
  
  m_objectives <- t(as.matrix(df_objectives))
  
  nondominated <- nondominated_points(m_objectives)
  
  rank <- nds_rank(m_objectives) #, partial)
  
  file_name_split <- unlist(strsplit(file_name, '_'))
  
  # Prepare plot for Pareto fronts 
  options(repr.plot.width=3, repr.plot.height=3)
  plot(c(), c(), 
       xlim=c(min(df_objectives$probability_complement), max(df_objectives$probability_complement)),
       ylim=c(min(df_objectives$distance_neg), max(df_objectives$distance_neg)),
       main=file_name_split[4],
       xlab='Probability Complement',
       ylab='Negated Pairwise Alignment Score')
  
  # Plot Pareto fronts
  for (idx in 1:max(rank)) {
    color <- 'blue'
    shape <- 1
    if (!idx %% 2) {
      color <- 'green'
      shape <- 4
    }
    
    # Prepare Pareto points
    pareto_front_points_indicies <- which(rank %in% c(idx))
    df_pareto_front_points <- df_objectives[pareto_front_points_indicies, 
                                            c('probability_complement', 'distance_neg')]
      
    # plots Pareto fronts 
    points(df_pareto_front_points$probability_complement, 
           df_pareto_front_points$distance_neg, 
           pch=shape)
    lines(df_pareto_front_points$probability_complement, 
          df_pareto_front_points$distance_neg, 
          col=color, lty=1)
  }
}

dev.off()
