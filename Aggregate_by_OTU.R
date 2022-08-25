library(argparse)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

################################
# Process command line arguments
################################

parser <- ArgumentParser(description="Plot mapped basepairs vs template coverage, optionally label reference community members (if present)")
parser$add_argument('-i', metavar='--input', type='character', 
                    help="Input file: input.tsv should be a tsv file generated with Compare_alignments.py")
parser$add_argument('-r', metavar='--reference_community', type='character', default="",
                    help=" Filepath to a reference community ;-separated csv, containing: Organism; Perc_gDNA")
parser$add_argument('-t', metavar='--taxonomic_level', type='character', default="species",
                    help="taxonomic level at which to compare alignments. One of [species, genus]")
parser$add_argument('-m', metavar='--experiment_metadata', type='character', default="",
                    help="experiment_metadata file: ;-separated csv file which should at least have two columns: Experiments; Alignment_database  Experiment and database should be named exactly as the output filenames")                    
args <- parser$parse_args()

filename <- args$i
taxlevel_species <- args$t=="species"
taxlevel_genus <- args$t=="genus"
reference <- args$r
metadata <- args$m


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
group_by(Experiment) %>%
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
            depth = bpTotal/total_template_length,
            p_readCount = sum(readCount, na.rm = TRUE)/unique(Total_readCount),
            readCount = sum(readCount, na.rm=TRUE),
            mean_readlength = bpTotal/readCount,
            refConsensusSum = sum(refConsensusSum),
            .groups = "drop")

if(metadata!=""){
  KMA <- merge(KMA, experiments, by="Experiment", all = TRUE)
}

if(reference!=""){
  KMA <- merge(KMA, ref_community, by.x="OTU", by.y = "Organism", all = TRUE)
}

output_file <- str_replace_all(filename, "\\.tsv", "_agg.tsv")
write.table(KMA, output_file, row.names=FALSE, sep="\t")