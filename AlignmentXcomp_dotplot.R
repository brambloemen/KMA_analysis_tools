.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)
library(data.table)


#the variable args below captures the arguments you pass from the 
#command line i.e. the names of the two files and stores them in a vector
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
##########################################
# functions to clean ARG and organism name
##########################################

clean_org_name <- function(organism){
  organism <- str_remove_all(organism, "Synthetic|\\[|\\]|\\scf.")
  organism <- str_replace_all(organism, "_", " ")
  organism <- str_replace_all(organism, "E.coli.*", "Escherichia coli")
  organism <- str_replace_all(organism, "veillonella", "Veillonella")
  organism <- str_replace_all(organism, "\\.", " ")
  organism <- str_extract(organism, "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?")
  organism <- str_replace(organism, "Clostridium difficil+e", "Clostridioides difficile")
  return(organism)
}

clean_ARG_names <- function(amr_gene){
  amr_gene <- str_remove(amr_gene, "(_|-).+")
  return(amr_gene)
}


######################################
# read in alignment comparison files
######################################

data <- fread(filename, sep = "\t", integer64 = "numeric")
data <- mutate(data, Method=str_remove(filename, ".tsv")) %>% as.data.frame()


##########################################
# Direct comparison between two alignments
##########################################

data <- mutate(data,
               Organism1 = clean_org_name(Template1),
               Organism1 = ifelse(is.na(Organism1), "Unmapped", Organism1),
               Organism2 = clean_org_name(Template2),
               Organism2 = ifelse(is.na(Organism2), "Unmapped", Organism2),
               Agree = Organism1 == Organism2) %>%
  group_by(Organism1, Organism2, Agree) %>%
  select(names(data[sapply(data, is.numeric)]))
data[is.na(data)] <- 0
data <- data %>% summarize_all(sum)
data <- mutate(data,
               Mean_mapped_bp = (n_match_bases1 + n_match_bases2)/2)

plot <- ggplot(data, aes(x = Organism1, y= Organism2)) +
  geom_point(aes(col=Agree, size=Mean_mapped_bp)) +
  theme_classic()


exportpath <- str_replace(filename, "\\.tsv?", "_Xdotplot.png")
print(exportpath)
png(exportpath)
print(plot)
dev.off()
