library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(optparse)

# Process input arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input file: summary data of KMA results, generated with Summarize_KMA_results.py", metavar="character"),
  make_option(c("-t", "--taxonomic_level"), type="character", default="species", 
              help="Taxonomic level at which to group [default= %default , genus]", metavar="character"),
  make_option(c("-r", "--reference_community"), type="character", default=NULL, 
              help="reference_community file: ;-separated csv file with two columns: Organism, Perc_gDNA", metavar="character"),
  make_option(c("-m", "--experiment_metadata"), type="character", default=NULL, 
              help="experiment_metadata file: ;-separated csv file which should at least have two columns:\n - Experiments\n - Alignment_database\n
              Experiment and database should be named exactly as the output filenames", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_file <- opt$input
taxlevel_species <- opt$taxonomic_level=="species"
taxlevel_genus <- opt$taxonomic_level=="genus"
ref_com_fp <- opt$reference_community=="Mean_mapped_bp"
experiments_fp <- opt$experiment_metadata=="Diff_mapped_bp"


######################
# Process input data 
######################

# mock community (reference)
ref_community <- read.csv(ref_com_fp, sep = ";", dec = ",")
ref_community_name <- str_extract(ref_com_fp, "/\\w+.csv")
ref_community_name <- str_remove_all(ref_community_name, "/|.csv")
ref_community <- ref_community %>%
  arrange(desc(Perc_gDNA), Organism)

# sort by descending zymo GMS gDNA% -> make organisms into ordered factors
ref_community$Organism <- factor(ref_community$Organism, levels = ref_community$Organism,
                           ordered = TRUE)

# Experiment metadata
experiments <- read.csv(experiments_fp, sep=";")
experiments$Experiment <- paste0(experiments$Experiment , "_", experiments$Alignment_database)

# KMA results
KMA <- fread(input_file, sep="\t", integer64 = "numeric")

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
KMA <- merge(KMA, experiments, by="Experiment", all = TRUE)
KMA <- merge(KMA, ref_community, by.x="OTU", by.y = "Organism", all = TRUE)

write.csv(KMA, "./")