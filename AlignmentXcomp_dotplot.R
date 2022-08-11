library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

# Process input arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input file comparing two alignments (Comparison_alignments.py)", metavar="character"),
  make_option(c("-t", "--taxonomic_level"), type="character", default="species", 
              help="Taxonomic level at which to group [default= %default , genus]", metavar="character"),
  make_option(c("-s", "--statistic"), type="character", default="Mean_mapped_bp", 
              help="Statistic to show on plot [default= %default , Diff_mapped_bp]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

filename <- opt$input
taxlevel_species <- opt$taxonomic_level=="species"
taxlevel_genus <- opt$taxonomic_level=="genus"
stat_mean <- opt$statistic=="Mean_mapped_bp"
stat_diff <- opt$statistic=="Diff_mapped_bp"


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
               Genus1 = str_extract(Organism1, "[:upper:]{1}[:lower:]+"),
               Genus2 = str_extract(Organism2, "[:upper:]{1}[:lower:]+"),
               Agree = Organism1 == Organism2)


if(taxlevel_species){
  data <- data %>%
    group_by(Organism1, Organism2, Agree) %>%
    select(names(data[sapply(data, is.numeric)]))
  data[is.na(data)] <- 0
  data <- data %>% summarize_all(sum)
  data <- mutate(data,
                 Mean_mapped_bp = (n_match_bases1 + n_match_bases2)/2,
                 Diff_mapped_bp = abs(n_match_bases1 - n_match_bases2),
                 Rel_diff_mapped_bp = Diff_mapped_bp/Mean_mapped_bp)

  if(stat_mean){
    plot <- ggplot(data, aes(x = Organism1, y= Organism2)) +
      geom_point(aes(col=Agree, size=Mean_mapped_bp)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      guides(color = guide_legend(override.aes = list(size = 8)))
  } else if(stat_diff){
    plot <- ggplot(data, aes(x = Organism1, y= Organism2)) +
      geom_point(aes(col=Agree, size=Rel_diff_mapped_bp)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      guides(color = guide_legend(override.aes = list(size = 8)))
  }
  
} else if(taxlevel_genus){
  data <- data %>%
    group_by(Genus1, Genus2, Agree) %>%
    select(names(data[sapply(data, is.numeric)]))
  data[is.na(data)] <- 0
  data <- data %>% summarize_all(sum)
  data <- mutate(data,
                 Mean_mapped_bp = (n_match_bases1 + n_match_bases2)/2,
                 Diff_mapped_bp = abs(n_match_bases1 - n_match_bases2),
                 Rel_diff_mapped_bp = Diff_mapped_bp/Mean_mapped_bp)
  
  if(stat_mean){
    plot <- ggplot(data, aes(x = Genus1, y= Genus2)) +
      geom_point(aes(col=Agree, size=Mean_mapped_bp)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      guides(color = guide_legend(override.aes = list(size = 8)))
  } else if(stat_diff){
    plot <- ggplot(data, aes(x = Genus1, y= Genus2)) +
      geom_point(aes(col=Agree, size=Rel_diff_mapped_bp)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      guides(color = guide_legend(override.aes = list(size = 8)))
  }
}

pdf(NULL)
exportpath <- str_replace(filename, "\\.tsv?", "_Xdotplot.png")
print(exportpath)
print(plot)
ggsave(exportpath, width = 16, height = 9, dpi = 100)

