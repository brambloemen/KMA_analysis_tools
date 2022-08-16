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
  reference <- ""
  metadata <- ""
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
      
    }else if (args[i] == "-r"){
      if (!str_detect(args[[i+1]], "\\.csv")){
        usage()
      }
      reference <- args[[i+1]]

    }else if (args[i] == "-m"){
      if (!str_detect(args[[i+1]], "\\.csv")){
        usage()
      }
      metadata <- args[[i+1]]
    }
    
  }
  argslist <- list(file=filename, taxlevel=taxlevel, ref=reference, metadata=metadata)
  return(argslist)
  
}

usage <- function() {
  cat(
    "usage: Rscript AlignmentXcomp_dotplot.R -i input.tsv -s statistic -t taxonomic_level",
    "  -i: Input file: tsv file with summary data of KMA results, generated with Summarize_KMA_results.py", 
    "  -r: Filepath to a reference ;-separated csv, containing 2 columns: Organism;Perc_gDNA, default: None",
    "  -t: taxonomic level at which to compare alignments. One of [species, genus]. Default: species",
    "  -m: experiment_metadata file: ;-separated csv file which should at least have two columns:\n - Experiments\n - Alignment_database\n
              Experiment and database should be named exactly as the output filenames", sep = "\n")
}

args <- main()
filename <- args$file
taxlevel_species <- args$taxlevel=="species"
taxlevel_genus <- args$taxlevel=="genus"
reference <- args$ref
metadata <- args$metadata


######################
# Process input data 
######################

# mock community (reference)
if(reference!=""){
  ref_community <- read.csv(reference, sep = ";", dec = ",")
  ref_community_name <- str_extract(reference, "/\\w+.csv")
  ref_community_name <- str_remove_all(ref_community_name, "/|.csv")
  ref_community <- ref_community %>%
    arrange(desc(Perc_gDNA), Organism)

  # sort by descending zymo GMS gDNA% -> make organisms into ordered factors
  ref_community$Organism <- factor(ref_community$Organism, levels = ref_community$Organism,
                           ordered = TRUE)
}

# Experiment metadata
if(metadata!=""){
  experiments <- read.csv(metadata, sep=";")
  experiments$Experiment <- paste0(experiments$Experiment , "_", experiments$Alignment_database)
}

# KMA results
KMA <- fread(filename, sep="\t", integer64 = "numeric")

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

KMA$OTU <- clean_org_name(KMA$`# refSequence`)
KMA <- KMA %>% 
  mutate(Total_bp = sum(bpTotal, na.rm = TRUE),
         Total_readCount = sum(readCount, na.rm = TRUE)) %>%
  group_by(OTU, Experiment) %>%
  summarize(p_bpTotal = sum(bpTotal, na.rm = TRUE)/unique(Total_bp),
            mean_queryID = weighted.mean(Query_Identity, bpTotal, na.rm=TRUE),
            mean_query_coverage = weighted.mean(Query_Coverage, bpTotal, na.rm=TRUE),
            mean_templateID = weighted.mean(Template_Identity, Template_length, na.rm=TRUE),
            total_template_length = sum(Template_length, na.rm=TRUE),
            mean_template_length = mean(Template_length, na.rm=TRUE),
            mean_template_coverage = weighted.mean(Template_Coverage, bpTotal, na.rm=TRUE),
            bpTotal = sum(bpTotal, na.rm = TRUE),
            p_readCount = sum(readCount, na.rm = TRUE)/unique(Total_readCount),
            readCount = sum(readCount, na.rm=TRUE),
            mean_readlength = bpTotal/readCount,
            refConsensusSum = sum(refConsensusSum),
            .groups = "drop")

if(metadata!="" & reference!=""){
  KMA <- merge(KMA, experiments, by="Experiment", all = TRUE)
  KMA <- merge(KMA, ref_community, by.x="OTU", by.y = "Organism", all = TRUE)
}

output_file <- str_replace_all(filename, "\\.tsv", "_agg.tsv")
write.csv(KMA, output_file, row.names=FALSE)