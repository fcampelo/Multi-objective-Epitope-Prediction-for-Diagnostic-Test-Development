# this script uses the devel-next version of the epitopes package:
# > devtools::install_github("fcampelo/epitopes@devel-next", force = TRUE)

library(epitopes)
library(dplyr)

to_gen <- data.frame(folder = c("01_EBV", "02_HepC", "03_Ovolvulus"),
                     taxids = c("10376", "11102", "6282"))

taxonomy <- readRDS("./taxonomy.rds")
peptides <- readRDS("./peptides.rds") %>%
  filter(Class == 1)

for (i in seq_along(to_gen$folder)){
  x <- peptides %>% 
    epitopes::filter_epitopes(removeIDs     = to_gen$taxids[i], 
                              orgID_column  = "Info_organism_id",
                              hostID_column = "Info_host_id", 
                              tax_list = taxonomy)
  saveRDS(x, paste0(to_gen$folder[i], "/epitopes_except_organism.rds"))
}
