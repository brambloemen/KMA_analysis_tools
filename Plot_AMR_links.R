library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(argparse)


################################
# Process command line arguments
################################

parser <- ArgumentParser(description="Plot mapped basepairs vs template coverage, optionally label reference community members (if present)")
parser$add_argument('-i', metavar='--input', type='character', 
                    help="Input file: input.tsv should be a tsv file generated with Compare_alignments.py")
parser$add_argument('-r', metavar='--reference_ARGprofile', type='character', default="",
                    help=" Filepath to a reference ;-separated csv, containing 2 columns: Organism;ARG, default: None")
parser$add_argument('-t', metavar='--taxonomic_level', type='character', default="species",
                    help="taxonomic level at which to compare alignments. One of [species, genus]")
args <- parser$parse_args()

filename <- args$i
taxlevel_species <- args$t=="species"
taxlevel_genus <- args$t=="genus"
reference <- args$r


###################################
# read Reference community AMR data
###################################

if(reference!=""){
  reference_data <- read.csv(reference, sep = ";") %>%
    mutate(Reference_AMR_link=TRUE,
           Species = Organism) %>%
    select(-Organism)
  reference_ARGs <- data.frame(reference_data$ARG) %>% mutate(Reference_AMR=TRUE)
}


##########################################
# functions to clean and summarize data  #
##########################################

clean_org_name <- function(organism){
  if(organism=="Unmapped"){
    organism <- "Unmapped"
  }
  else{
    organism <- str_remove_all(organism, "Synthetic|\\[|\\]|\\scf.")
    organism <- str_replace_all(organism, "_", " ")
    organism <- str_replace_all(organism, "E.coli.*", "Escherichia coli")
    organism <- str_replace_all(organism, "veillonella", "Veillonella")
    organism <- str_replace_all(organism, "\\.", " ")
    organism <- str_extract(organism, "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?")
    organism <- str_replace(organism, "Clostridium difficil+e", "Clostridioides difficile")
  }
  
  return(organism)
}

clean_ARG_names <- function(amr_gene){
  amr_gene <- str_remove(amr_gene, "_.+")
  amr_gene <- str_remove(amr_gene, "-\\d+.+$")
  return(amr_gene)
}

aggregate_by_OTU <- function(data){
  data <- data %>%
    summarize(n_reads=sum(n_reads, na.rm = TRUE),
              Total_readlength=sum(Total_readlength, na.rm = TRUE),
              n_match_bases1=sum(n_match_bases1, na.rm = TRUE),
              n_match_bases2=sum(n_match_bases2, na.rm = TRUE),
              Template1_length=sum(Template1_length)/n(),
              Template2_length=sum(Template2_length)/n(),
              Organism_QID = n_match_bases1/Total_readlength,
              ARG_TID = n_match_bases2/(Template2_length*n_reads))
    return(data)
}


#######################
# plotting functions  #
#######################

plot_links_byreads <- function(data, taxa){
  plot <- data %>%
      ggplot(aes(x=ARG, y={{taxa}})) +
      geom_point(aes(size=AMR_links)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      guides(color = guide_legend(override.aes = list(size = 8)))
  return(plot)
} 


####################################
# read in alignment comparison file
####################################

data <- fread(filename, sep = "\t", integer64 = "numeric", fill=TRUE)

data <- mutate(data, 
               Species=sapply(Template1, clean_org_name),
               Genus=str_extract(Species, "([:upper:]{1}[:lower:]+)|Unmapped"), 
               ARG=clean_ARG_names(Template2)) %>%
    filter(ARG != "") %>%
    filter(!is.na(ARG) & ARG != "Unmapped") 

if(taxlevel_species){
  
  data <- data %>%
    group_by(Species, ARG) %>%
    aggregate_by_OTU() %>%
    mutate(AMR_links = case_when(
                Species=="Unmapped" ~  Total_readlength * ARG_TID,
                TRUE ~ Total_readlength * Organism_QID * ARG_TID)) %>%
    filter(AMR_links > 0)
  
  if(reference==""){
    plot <- plot_links_byreads(data, Species)
  }else{
    data <- merge(data, reference_data, by=c("Species", "ARG"), all=TRUE) %>%
      mutate(Reference_AMR_link=ifelse(Reference_AMR_link, Reference_AMR_link, FALSE),
             Reference_AMR_link=ifelse(is.na(Reference_AMR_link), FALSE, Reference_AMR_link)) 
    
    plot <- plot_links_byreads(data, Species) + geom_point(aes(col=Reference_AMR_link, size=AMR_links))
  }
  
} else if(taxlevel_genus){
  
  data <- data %>%
    group_by(Genus, ARG) %>%
    aggregate_by_OTU() %>%
    mutate(AMR_links = case_when(
                Genus=="Unmapped" ~  Total_readlength * ARG_TID,
                TRUE ~ Total_readlength * Organism_QID * ARG_TID)) %>%
    filter(AMR_links > 0)
  
  if(reference==""){
    plot <- plot_links_byreads(data, Genus)
  }else{
    reference_data <- mutate(reference_data, Genus=str_extract(Species, "[:upper:]{1}[:lower:]+"))
    
    data <- merge(data, reference_data, by=c("Genus", "ARG"), all=TRUE) %>%
      mutate(Reference_AMR_link=ifelse(Reference_AMR_link, Reference_AMR_link, FALSE),
             Reference_AMR_link=ifelse(is.na(Reference_AMR_link), FALSE, Reference_AMR_link)) 
    
    plot <- plot_links_byreads(data, Genus) + geom_point(aes(col=Reference_AMR_link, size=AMR_links))
  }
}

pdf(NULL)
exportpath <- str_replace(filename, "\\.tsv?", "_AMRlinks.png")
print(exportpath)
print(plot)
ggsave(exportpath, width = 16, height = 9, dpi = 100)