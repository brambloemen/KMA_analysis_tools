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
parser$add_argument('-s', metavar='--statistic', type='character', default="Mean_mapped_bp",
                    help="statistic to plot, one of: [Mean_mapped_bp: mean of mapped bp to a given template of the two compared alignments, Diff_mapped_bp: difference in mapped bp] Default:Mean_mapped_bp")
parser$add_argument('-t', metavar='--taxonomic_level', type='character', default="species",
                    help="taxonomic level at which to compare alignments. One of [species, genus]")
args <- parser$parse_args()

filename <- args$i
taxlevel_species <- args$t=="species"
taxlevel_genus <- args$t=="genus"
stat_mean <- args$s=="Mean_mapped_bp"
stat_diff <- args$s=="Diff_mapped_bp"


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
                 Mean_mapped_bp = as.numeric(n_match_bases1 + n_match_bases2)/2,
                 Diff_mapped_bp =  as.numeric(abs(n_match_bases1 - n_match_bases2)),
                 Rel_diff_mapped_bp = Diff_mapped_bp/Mean_mapped_bp)

  if(stat_mean){
    plot <- ggplot(data, aes(x = Organism1, y= Organism2)) +
      geom_point(aes(col=Agree, size=Mean_mapped_bp)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      labs(x= "Alignment 1", y= "Alignment 2") +
      guides(color = guide_legend(override.aes = list(size = 8)))
  } else if(stat_diff){
    plot <- ggplot(data, aes(x = Organism1, y= Organism2)) +
      geom_point(aes(col=Agree, size=Rel_diff_mapped_bp)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      labs(x= "Alignment 1", y= "Alignment 2") +
      guides(color = guide_legend(override.aes = list(size = 8)))
  }
  
} else if(taxlevel_genus){
  data <- data %>%
    group_by(Genus1, Genus2, Agree) %>%
    select(names(data[sapply(data, is.numeric)]))
  data[is.na(data)] <- 0
  data <- data %>% summarize_all(sum)
  data <- mutate(data,
                 Mean_mapped_bp =  as.numeric(n_match_bases1 + n_match_bases2)/2,
                 Diff_mapped_bp =  as.numeric(abs(n_match_bases1 - n_match_bases2)),
                 Rel_diff_mapped_bp = Diff_mapped_bp/Mean_mapped_bp)
  
  if(stat_mean){
    plot <- ggplot(data, aes(x = Genus1, y= Genus2)) +
      geom_point(aes(col=Agree, size=Mean_mapped_bp)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      labs(x= "Alignment 1", y= "Alignment 2") +
      guides(color = guide_legend(override.aes = list(size = 8)))
  } else if(stat_diff){
    plot <- ggplot(data, aes(x = Genus1, y= Genus2)) +
      geom_point(aes(col=Agree, size=Rel_diff_mapped_bp)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90)) +
      labs(x= "Alignment 1", y= "Alignment 2") +
      guides(color = guide_legend(override.aes = list(size = 8)))
  }
}

pdf(NULL)
exportpath <- str_replace(filename, "\\.tsv?", "_Xdotplot.png")
print(exportpath)
print(plot)
ggsave(exportpath, width = 16, height = 9, dpi = 100)

