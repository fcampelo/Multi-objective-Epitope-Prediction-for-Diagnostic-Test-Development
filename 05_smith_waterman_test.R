library(BiocManager)
library(Biostrings)


#s1 <- AAString("FSNRDVYDSGSGVNNYDGDY")
#s2 <- AAString("SFSNRDVYDSGSGVNNYDGDY")

#s1 <- AAString("GMTPLTEEEYQYSLQW")
#s2 <- AAString("MTPLTEEEYQYSLQW")

#s1 <- AAString("PTTGETAE")
#s2 <- AAString("KPTTEGGVETE")

#s1 <- AAString("MINS")
#s2 <- AAString("MAKLILVS")

#s1 <- AAString("DLTDIKRN")
#s2 <- AAString("DLNTLDINE")

s1 <- AAString("ERSVSQCKTSDSRNRSSSRRNSPEYDISQITEQESEKS")
s2 <- AAString("PKYVKQNTLKLATGMRNVPEKQT") #EYEEDRGKGFFPQTITQEYE")

#s1 <- AAString("SKSRSESLVVSRQ")
#s2 <- AAString("GKGRGLSLSR")

#s1 <- AAString("RSEKSSVREF")
#s2 <- AAString("RTSITNTGLT")




## Align two amino acid sequences with the BLOSUM62 matrix

#pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM62", gapOpening = 3, gapExtension = 1)

## See how the gap penalty influences the alignment
#pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM62", gapOpening = 6, gapExtension = 2)

## See how the substitution matrix influences the alignment
#pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM50", gapOpening = 3, gapExtension = 1)

#pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM45", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM50", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM62", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM80", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM100", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "PAM30", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "PAM40", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "PAM70", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "PAM120", gapOpening = 3, gapExtension = 1)
#pairwiseAlignment(s1, s2, substitutionMatrix = "PAM250", gapOpening = 3, gapExtension = 1)

pairwiseAlignment(s1, s2, type = "overlap", substitutionMatrix = "BLOSUM45", gapOpening = 3, gapExtension = 1)
pairwiseAlignment(s1, s2, type = "overlap", substitutionMatrix = "PAM30", gapOpening = 3, gapExtension = 1)



#if (interactive()) {
#  ## Compare our BLOSUM62 with BLOSUM62 from ftp://ftp.ncbi.nih.gov/blast/matrices/
#  data(BLOSUM62)
#  BLOSUM62["Q", "Z"]
#  file <- "ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62"
#  b62 <- as.matrix(read.table(file, check.names=FALSE))
#  b62["Q", "Z"]
#}


#data(PAM30)
#data(PAM250)

#data(BLOSUM45)
