require("dplyr")
require("emoa")
require("mco")
require("ggplot2")
library("plotly")

PLOT_SEQUENCE <- FALSE
file_name <- './data/predictions_pairwiseAlignment_PAM30_local_2021-07-05_16_47_49.csv'

df <- read.csv(file=file_name)

df_complemented_negated <- df %>%
  mutate(probability_complement = 1 - Probability) %>%
  select(Protein, Sequence, 
         Probability, pair_align_score, 
         probability_complement, 
         Info_peptide, Info_PepID) #%>% head() 

df_objectives <- df_complemented_negated %>%
  select(probability_complement, pair_align_score) #%>% head() 

m_objectives <- t(as.matrix(df_objectives))

nondominated <- nondominated_points(m_objectives)

rank <- nds_rank(m_objectives) #, partial)

file_name_split <- unlist(strsplit(file_name, '_'))

if (PLOT_SEQUENCE) {
  pl <- ggplot(data=df_complemented_negated, aes(x=Probability, y=pair_align_score)) +
    ggtitle("PAM30") +
    geom_point(aes(text=sprintf("Protein ID: %s Sequence: %s<br>Info_PepID: %s Info_peptide: %s", 
                                Protein, Sequence, 
                                Info_PepID, Info_peptide)))
} else {
  pl <- ggplot(data=df_complemented_negated, aes(x=Probability, y=pair_align_score)) +
    ggtitle("PAM30") +
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

chart_path = paste('./output/chart_', file_name_split[5], file_name_split[6], file_name_split[7], sep='')

htmlwidgets::saveWidget(chart, selfcontained=TRUE, paste(chart_path, '.html', sep=''))
