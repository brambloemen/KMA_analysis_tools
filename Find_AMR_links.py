import pysam
import pandas as pd
import argparse
import re
import csv
from collections import defaultdict
import time
import logging
import sys

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')
parser = argparse.ArgumentParser(description='Retrieve alignments from alignment file 1 in alignment file 2; \nPrint results to standard output')
parser.add_argument("-t", metavar="taxa.bam", dest="bam_taxa", type=str, help="Filepath of bam with alignment to taxa database")
# parser.add_argument("-r", dest="taxa_res", type=str, help="Filepath of res file of KMA alignmnent to taxa database")
parser.add_argument("-a", metavar="AMR.bam", dest="AMR_bam", type=str, help="Filepath of bam with alignments to resistance database")
parser.add_argument("--tc", type=float, help="Minimum template coverage for a taxa template to be accepted", default=0.0)
parser.add_argument("--pid_arg", type=float, help="Minimum Percentage Identity to AMR gene template", default=0.0)
args = parser.parse_args()

# Filter for taxa
fp_taxa_res = re.sub("\..am", ".res", args.bam_taxa)
Filtered_taxa = pd.read_csv(fp_taxa_res, sep="\t")
Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_Coverage"] >= args.tc]["#Template"]
Filtered_taxa = Filtered_taxa.tolist()

# Filter for AMR genes
fp_AMR_res = re.sub("\..am", ".res", args.AMR_bam)
Filtered_AMR = pd.read_csv(fp_AMR_res, sep="\t")
Filtered_AMR = Filtered_AMR[Filtered_AMR["Template_Identity"] >= args.pid_arg]["#Template"]
Filtered_AMR = Filtered_AMR.tolist()


"""
Loop over alignment to taxa
"""
logging.info(f"Parsing alignment 1: {args.bam_taxa}")
align1 = pysam.AlignmentFile(args.bam_taxa, "rb", threads=1)
align1_refseq = {ref["SN"]:ref["LN"] for ref in align1.header.to_dict()["SQ"]}
align1_n_lines = 0

# align1 dictionary: {read:[Template1, readlength, matching bases, Template1 total reference length]}

align1_reads = defaultdict(lambda: [str, 0, 0, 0])
for read in align1:
    align1_n_lines += 1
    # skip alignment if total template coverage was too low
    if read.reference_name not in Filtered_taxa:
        continue
    # store read identifiers
    n = re.search("read=\d+\sch=\d+",read.query_name).group()
    n = re.sub("\D", "", n)
    Template1 = read.reference_name
    l = read.query_length
    cigarEQ = read.get_cigar_stats()[0][7]
    Template1_l = align1_refseq.get(Template1, None)

    align1_reads[n] = [Template1, l, cigarEQ, Template1_l]
    
logging.info(f"Done parsing alignment 1. Processed {align1_n_lines} alignments")


"""
Loop over alignment to resistance database
"""

# summary dictionary. structure: {AMR gene: {Species: [#reads, total readlength, #exactly matched bases, readlength mapped to AMR, exactly matched bases to AMR, AMR template length]}}
alignment_links = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0, 0, 0]))
logging.info(f"Parsing alignment 2: {args.AMR_bam}")
align2 = pysam.AlignmentFile(args.AMR_bam, "rb", threads=1)
align2_refseq = {ref["SN"]:ref["LN"] for ref in align2.header.to_dict()["SQ"]}
align2_n_lines = 0

for read in align2:
    align2_n_lines += 1
    # read identifiers
    n = re.search("read=\d+\sch=\d+",read.query_name).group()
    n = re.sub("\D", "", n)
    # alignment stats
    if read.reference_name == None: # not aligned to AMR --> skip
        continue
    if read.reference_name not in Filtered_AMR: # KMA alignment to AMR gene template was not identical enough
        continue
    else:
        Template2 = read.reference_name
    l = read.query_length
    cigarEQ = read.get_cigar_stats()[0][7] # EQ: exact sequence match (base in query matches base in template)
    Template2_l = align2_refseq.get(Template2, None)
    # matching bases/AMR length == Percentage Identity --> should also be filtered, !! currently does not consider reads that align to only part of the AMR
    if 100*cigarEQ/Template2_l < args.pid_arg:
        continue


    if n in align1_reads:
        Template1, Template1_readl, Template1_cigarEQ, Template1_l = align1_reads.pop(n) # pop matched reads from amr_reads: unmapped reads will remain
        alignment_links[Template1][Template2][0] += 1
        alignment_links[Template1][Template2][1] += Template1_readl # will be total length of all reads aligned to both template 1 and 2
        alignment_links[Template1][Template2][2] += Template1_cigarEQ
        alignment_links[Template1][Template2][3] += cigarEQ 
        alignment_links[Template1][Template2][4] = Template1_l
        alignment_links[Template1][Template2][5] = Template2_l
    else:
        alignment_links["Unmapped"][Template2][0] += 1
        alignment_links["Unmapped"][Template2][1] += l # will be total length of all reads aligned to template 2, with no match for temp
        alignment_links["Unmapped"][Template2][2] = None
        alignment_links["Unmapped"][Template2][3] += cigarEQ 
        alignment_links["Unmapped"][Template2][4] = None
        alignment_links["Unmapped"][Template2][5] = Template2_l

logging.info(f"Done parsing alignment 2. Processed {align2_n_lines} alignments")


# print results to standard output
logging.info(f"Writing results to tsv")
header = ["Template1", "Template2", "n_reads", "Total_readlength", "n_match_bases1", "n_match_bases2", "Template1_length", "Template2_length"]
writer = csv.writer(sys.stdout, delimiter="\t")
writer.writerow(header)

for temp1, temp2 in alignment_links.items():
    for template2, stats in temp2.items():
        writer.writerow([temp1, template2] + stats)

logging.info(f"Done")