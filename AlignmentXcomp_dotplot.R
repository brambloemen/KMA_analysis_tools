library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

################################
# Process command line arguments
################################
main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  #default settings
  stat <- "Mean_mapped_bp"
  taxlevel <- "species"
  
  for (i in 1:length(args)){
    
    if (args[i] == "-i"){
      if (!str_detect(args[[i+1]], "\\.tsv")){
        usage()
      }
      filename <- args[[i+1]]
    }
    else if (args[i] == "-t"){
      if (!(args[[i+1]] %in% c("species", "genus"))){
        usage()
      }
      if (args[[i+1]]=="genus"){
        taxlevel <- "genus"
      }else {
        taxlevel <- "species"
      }
      
    }else if (args[i] == "-s"){
      if (!(args[[i+1]] %in% c("Mean_mapped_bp", "Diff_mapped_bp", ""))){
        usage()
      }
      if (args[[i+1]]=="Diff_mapped_bp"){
        stat <- "Diff_mapped_bp"
      }
      else {
        stat <- "Mean_mapped_bp"
      }
    }
    
  }
  argslist <- list(file=filename, taxlevel=taxlevel, stat=stat)
  return(argslist)
  
}

usage <- function() {
  cat(
    "usage: Rscript AlignmentXcomp_dotplot.R -i input.tsv -s statistic -t taxonomic_level",
    "  -i: input.tsv should be a tsv file generated with Compare_alignments.py", 
    "  -s: statistic to plot, one of: \n[Mean_mapped_bp: mean of mapped bp to a given template of the two compared alignments, Diff_mapped_bp: difference in mapped bp]\n Default:Mean_mapped_bp",
    "  -t: taxonomic level at which to compare alignments. One of [species, genus]",sep = "\n")
}

args <- main()
filename <- args$file
taxlevel_species <- args$taxlevel=="species"
taxlevel_genus <- args$taxlevel=="genus"
stat_mean <- args$stat=="Mean_mapped_bp"
stat_diff <- args$stat=="Diff_mapped_bp"


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

