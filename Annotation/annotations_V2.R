library(openxlsx)
library(readxl)

# import annotations from different sources
circRNA <- read.xlsx("circRNA_annotation.xlsx")
head(circRNA)
plosOne <-  read_xls("journal.pone.0187581_annotation.xls")
head(plosOne)
frontiers <- read_xls("frontiers.2019.00526_annotation.xls")
head(frontiers)
scirep <- read.xlsx("scirep.47189_annotation.xlsx")
head(scirep)

# import complete list of probes
arrayStar <- read.delim("A-GEOD-21825_clean.txt")
head(arrayStar)

# create a basis for annotations file
annotations <- rbind(circRNA, plosOne, frontiers, scirep)

#names of data
data_names <- colnames(annotations)

# make a new dataframe to be able to compare it to original
arrayStar1 <- arrayStar

# adding annotations
data_names
arrayStar1$Alias <- annotations[match(arrayStar1$ProbeName, annotations$probeID),3]
arrayStar1$chrom <- annotations[match(arrayStar1$ProbeName, annotations$probeID),4]
arrayStar1$strand <- annotations[match(arrayStar1$ProbeName, annotations$probeID),5]
arrayStar1$txStart <- annotations[match(arrayStar1$ProbeName, annotations$probeID),6]
arrayStar1$txEnd <- annotations[match(arrayStar1$ProbeName, annotations$probeID),7]
arrayStar1$circRNA_type <- annotations[match(arrayStar1$ProbeName, annotations$probeID),8]
arrayStar1$best_transcript <- annotations[match(arrayStar1$ProbeName, annotations$probeID),9]
arrayStar1$GeneSymbol <- annotations[match(arrayStar1$ProbeName, annotations$probeID),10]
arrayStar1$Sequence <- annotations[match(arrayStar1$ProbeName, annotations$probeID),11]

# check imported data
head(arrayStar1)
sum(is.na(arrayStar1$chrom)) # how many probes do not have annotations

# export probes & annotations
write.table(arrayStar1, file = "complete_annotations_V2.txt", sep='\t', row.names = FALSE)
