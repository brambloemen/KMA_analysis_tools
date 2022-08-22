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
parser$add_argument('-t', metavar='--taxonomic_level', type='character', default="Mean_mapped_bp",
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
  reference_data <- read.csv(reference, sep = "\t") %>%
    mutate(Reference_AMR_link=TRUE)
  reference_ARGs <- data.frame(reference_data$ARG) %>% mutate(Reference_AMR=TRUE)
}


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


####################################
# read in alignment comparison file
####################################

data <- fread(filename, sep = "\t", integer64 = "numeric")

data <- mutate(data, 
               Template1=clean_org_name(Template1),
               Organism=Template1,
               Genus=str_extract(Organism, "[:upper:]{1}[:lower:]+"),
               Template2=clean_ARG_names(Template2),
               ARG=Template2)

if(taxlevel_species){
  
  data <- data %>%
    group_by(Organism, ARG) %>%
    summarize(n_reads=sum(n_reads, na.rm = TRUE),
              sum_readlength=sum(sum_readlength, na.rm = TRUE),
              n_match_bases1=sum(n_match_bases1, na.rm = TRUE),
              n_match_bases2=sum(n_match_bases2, na.rm = TRUE),
              Template1_length=sum(Template1_length)/n(),
              Template2_length=sum(Template2_length)/n()) %>%
    mutate(Coverage_ARG = n_match_bases2)
  
  if(reference==""){
    plot <- data %>%
      ggplot(aes(x=ARG, y=Organism)) +
      geom_point(aes(size=n_reads)) +
      theme_classic() +
      guides(color = guide_legend(override.aes = list(size = 8)))
  }else{
    data <- merge(data, reference_data, by=c("Organism", "ARG"), all=TRUE) %>%
      mutate(Reference_AMR_link=ifelse(Reference_AMR_link, Reference_AMR_link, FALSE),
             Reference_AMR_link=ifelse(is.na(Reference_AMR_link), FALSE, Reference_AMR_link)) 
    
    plot <- data %>%
      ggplot(aes(x=ARG, y=Organism)) +
      geom_point(aes(col=Reference_AMR_link, size=n_reads)) +
      theme_classic() +
      guides(color = guide_legend(override.aes = list(size = 8)))
  }
  
} else if(taxlevel_genus){
  
  data <- data %>% 
    group_by(Genus, ARG) %>%
    summarize(n_reads=sum(n_reads, na.rm = TRUE),
              sum_readlength=sum(sum_readlength, na.rm = TRUE),
              n_match_bases1=sum(n_match_bases1, na.rm = TRUE),
              n_match_bases2=sum(n_match_bases2, na.rm = TRUE),
              Template1_length=sum(Template1_length)/n(),
              Template2_length=sum(Template2_length)/n()) %>%
    mutate(Coverage_ARG = n_match_bases2)
  
  if(reference==""){
    plot <- data %>%
      ggplot(aes(x=ARG, y=Genus)) +
      geom_point(aes(size=n_reads)) +
      theme_classic() +
      guides(color = guide_legend(override.aes = list(size = 8)))
  }else{
    reference_data <- mutate(reference_data, Genus=str_extract(Organism, "[:upper:]{1}[:lower:]+"))
    
    data <- merge(data, reference_data, by=c("Genus", "ARG"), all=TRUE) %>%
      mutate(Reference_AMR_link=ifelse(Reference_AMR_link, Reference_AMR_link, FALSE),
             Reference_AMR_link=ifelse(is.na(Reference_AMR_link), FALSE, Reference_AMR_link)) 
    
    plot <- data %>%
      ggplot(aes(x=ARG, y=Genus)) +
      geom_point(aes(col=Reference_AMR_link, size=n_reads)) +
      theme_classic() +
      guides(color = guide_legend(override.aes = list(size = 8)))
  }
}

pdf(NULL)
exportpath <- str_replace(filename, "\\.tsv?", "_AMRlinks.png")
print(exportpath)
print(plot)
ggsave(exportpath, width = 16, height = 9, dpi = 100)