
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(argparse)
library(data.table)

################################
# Process command line arguments
################################
parser <- ArgumentParser(description="Plot expected versus observed relative abundance (mapped basepairs)")
parser$add_argument('-i', metavar='--input', type='character', 
                    help="Input file: KMA summary .tsv file, generated by Aggregate_by_OTU.R, with reference community data")
parser$add_argument('-c', metavar='--category', type='character', default="Experiment", 
                    help="Categorical variable to compare")

args <- parser$parse_args()

KMA <- read.csv(args$i, sep = "\t")
Compare_category <- args$c

#Graph theoretical comparison
graph <- function(df, categorical=!!sym(Compare_category)){
  ggplot(df, aes(fill=OTU, x={{categorical}}, y=p_bpTotal)) +
    geom_bar(position="fill", stat="identity", width = 0.9) +
    scale_x_discrete(expand = c(0, 0.6)) +
    scale_y_continuous(expand = c(0, 0),labels = scales::label_percent(accuracy = 1)) +
    ylab("Microbial Composition") + xlab("") +
    theme(axis.text.x = element_text(angle=-45, hjust=0.0, size=16),
          axis.text.y = element_text(size=16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          legend.position = "right",
          legend.title = element_text(size=20),
          legend.text = element_text(size=16),
          plot.margin = margin(2, 5.5, 2, 2, "cm")) +
    guides(fill = guide_legend(ncol = 1))
}


KMA <- KMA %>% filter(!is.na(Perc_gDNA))
ref_com <- KMA %>% select(OTU, Perc_gDNA) %>% 
  distinct() %>%
  arrange(desc(Perc_gDNA))
ref_com$p_bpTotal <- ref_com$Perc_gDNA/100
ref_com[Compare_category] = "Theoretical distribution"
KMA <- rbindlist(list(KMA, ref_com), fill = TRUE)

ref_OTU <- ordered(ref_com$OTU)
KMA$OTU <- ordered(KMA$OTU, levels=ref_OTU)

colors_needed <- length(ref_OTU)
palette <- c("#F1975A", "#B5B5B5", "#FFCD33", "#7CAFDD", "#997300", "#255E91", "#43682B", "#698ED0",
             "#A4A4A4", "#ED7D31", "#616161", "#264478", "#70AD47", "#5B9BD5", "#9E480E", "#E1C200",
             "#4170C4","#BF0000")
colors <- replicate(ceiling(colors_needed/(length(palette))), palette)
plot <- graph(KMA) + scale_fill_manual(values=colors)

pdf(NULL)
exportpath <- str_replace(args$i, "\\.tsv?", "_RA_v_refcom.png")
print(exportpath)
print(plot)
ggsave(exportpath, width = 16, height = 16, dpi = 100)

