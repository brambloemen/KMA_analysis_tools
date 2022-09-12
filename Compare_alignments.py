import pysam
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
parser.add_argument("bam1", metavar="--input_1.bam", type=str, help="Filepath of bam alignment 1")
parser.add_argument("bam2", metavar="--input_2.bam", type=str, help="Filepath of bam alignment 2")
parser.add_argument("--threads", type=int, help="Number of threads used to process SAM/BAM file", default=1)
args = parser.parse_args()


"""
Loop over alignment 1
"""
logging.info(f"Parsing alignment 1: {args.bam1}")
align1 = pysam.AlignmentFile(args.bam1, "rb", threads=args.threads)
align1_refseq = {ref["SN"]:ref["LN"] for ref in align1.header.to_dict()["SQ"]}
align1_n_lines = 0

# align1 dictionary: {read:[Template1, readlength, matching bases, Template1 total reference length]}

align1_reads = defaultdict(lambda: [str, 0, 0, 0])
for read in align1:
    # store read identifiers
    n = read.query_name
    Template1 = read.reference_name
    l = read.query_length
    cigarEQ = read.get_cigar_stats()[0][7]
    Template1_l = align1_refseq.get(Template1, None)

    align1_reads[n] = [Template1, l, cigarEQ, Template1_l]
    align1_n_lines += 1
logging.info(f"Done parsing alignment 1. Processed {align1_n_lines} alignments")


"""
Loop over alignment 2
"""

# summary dictionary. structure: {AMR gene: {Species: [#reads, total readlength, #exactly matched bases, readlength mapped to AMR, exactly matched bases to AMR, AMR template length]}}
alignment_links = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0, 0, 0]))
logging.info(f"Parsing alignment 2: {args.bam2}")
align2 = pysam.AlignmentFile(args.bam2, "rb", threads=args.threads)
align2_refseq = {ref["SN"]:ref["LN"] for ref in align2.header.to_dict()["SQ"]}
align2_n_lines = 0

for read in align2:
    # read identifiers
    n = read.query_name
    # alignment stats
    Template2 = read.reference_name
    l = read.query_length
    cigarEQ = read.get_cigar_stats()[0][7] # EQ: exact sequence match (base in query matches base in template)
    Template2_l = align2_refseq.get(Template2, None)


    if n in align1_reads:
        Template1, Template1_readl, Template1_cigarEQ, Template1_l = align1_reads.pop(n, None) # pop matched reads from amr_reads: unmapped reads will remain
        alignment_links[Template1][Template2][0] += 1
        alignment_links[Template1][Template2][1] += Template1_readl # will be total length of all reads aligned to both template 1 and 2
        alignment_links[Template1][Template2][2] += Template1_cigarEQ
        alignment_links[Template1][Template2][3] += cigarEQ 
        alignment_links[Template1][Template2][4] = Template1_l
        alignment_links[Template1][Template2][5] = Template2_l

    align2_n_lines += 1
logging.info(f"Done parsing alignment 2. Processed {align2_n_lines} alignments")

# print results to standard output
logging.info(f"Writing results to tsv")
header = ["Template1", "Template2", "n_reads", "Total_readlength", "n_match_bases1", "n_match_bases2", "Template1_length", "Template2_length"]
writer = csv.writer(sys.stdout, delimiter="\t")
writer.writerow(header)

for temp1, temp2 in alignment_links.items():
    for template2, stats in temp2.items():
        writer.writerow([temp1, template2] + stats)

# reads from bam 2 not recognized in bam 1
unmapped = defaultdict(lambda: {"Unmapped":[0, 0, 0, 0, 0, 0]})
# summarize remaining (= unmapped) alignments from align1_reads dictionary
for stats in align1_reads.values():
    Template1, Template1_readl, Template1_cigarEQ, Template1_l = stats
    unmapped[Template1]["Unmapped"][0] += 1
    unmapped[Template1]["Unmapped"][1] += Template1_readl # will be total length of all reads aligned to both template 1 and 2
    unmapped[Template1]["Unmapped"][2] += Template1_cigarEQ
    unmapped[Template1]["Unmapped"][3] = None
    unmapped[Template1]["Unmapped"][4] = Template1_l
    unmapped[Template1]["Unmapped"][5] = None
# write unmapped row per amr gene
for temp1, unmap_stats in unmapped.items():
    writer.writerow([temp1] + list(unmap_stats))

logging.info(f"Done")
