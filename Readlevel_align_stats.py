import pysam
import argparse
import re
import csv
from collections import defaultdict
import time
import logging

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')
parser = argparse.ArgumentParser(description='Extract read-level and template-level mapping statistics from KMA-analyes, including unmapped reads\n (which are often not returned by KMA).')
parser.add_argument("bam", metavar="--input.bam", type=str, help="Filepath of bamfile")
parser.add_argument("fastq", metavar="--input.fastq", type=str, help="Filepath of original fastq file on which KMA was performed")
"""
Note: In KMA versions > 1.4.3, the problem of reads missing from sam alignment files should be solved.
    --> This means that parsing the FastQ to retrieve missing reads is unnecessary.
    --> make the FastQ input optional
"""
args = parser.parse_args()
print(args)

t1=time.perf_counter()

readstatsfp = re.sub(".bam|.sam", "", args.bam) + "_readstats.tsv"
mapstatsfp = re.sub(".bam|.sam", "", args.bam) + "_mapstat.tsv"
with open(readstatsfp, "w") as csv_reads, open(mapstatsfp, "w") as csv_mapstat:

    csv_readstats = csv.writer(csv_reads, delimiter="\t")
    csv_readstats.writerow(["Read_ID", "Template", "Length", "Matched_bp"])

    csv_templatestats = csv.writer(csv_mapstat, delimiter="\t")
    csv_templatestats.writerow(["Template", "Readcount", "Total_readlength", "Total_bp"])

    # structure: [number of reads, combined length of reads, total sequence match bp]
    species_sum = defaultdict(lambda: [0, 0, 0])

    """
    First create a dict from the original fastq file: readID: readlength
    """
    logging.info("Parsing FASTQ")
    fastqfp = args.fastq
    fastq_reads = defaultdict(int)
    fqlines = 0

    with open(fastqfp) as fq:
        i = 0
        for line in fq:
            i += 1
            fqlines += 1
            if i == 1:
                seqid = re.search("read=\d+\sch=\d+", line).group() # smallest unique read id: read nr + channel number
                seqid = re.sub("\D", "", seqid)
            if i == 2:
                seqln = len(line.rstrip())
            if i == 4:
                fastq_reads[seqid] = seqln
                i = 0

    print(f"Total fastq reads: {fqlines/4}")
    t2=time.perf_counter()
    print(f"Processed {fqlines} lines in {t2 -t1:0.4f} seconds")

    """
    Next, loop over samfile. Write required alignment info to csv files.
    If sam read occurs in fastq, delete it from the fastq dictionary
    """
    logging.info("Parsing BAM file")
    samfile = pysam.AlignmentFile(args.bam, "rb", threads=1)
    samlines = 0
    for read in samfile:
        # Most simple Read ID
        n = re.search("read=\d+\sch=\d+",read.query_name).group()
        n = re.sub("\D", "", n)

        # remove sam read from fastq, default to None (don't throw keyerror)
        fastq_reads.pop(n, None)

        # alignment stats
        s = read.reference_name
        l = read.query_length
        # cigarM += read.get_cigar_stats()[0][0] # alignement match (a base in query aligns with a base in template)
        cigarEQ = read.get_cigar_stats()[0][7] # EQ: exact sequence match (base in query matches base in template)

        # summarize stats by species
        species_sum[s][0] += 1
        species_sum[s][1] += l
        species_sum[s][2] += cigarEQ

        csv_readstats.writerow([n, s, l, cigarEQ])
        samlines += 1

    # add Unmapped read stats from reads only present in fastq
    species_sum["Unmapped"][0] += len(fastq_reads) # number of reads
    species_sum["Unmapped"][1] += sum(fastq_reads.values()) # total read lengths = sum of lengths
    species_sum["Unmapped"][2] += sum(fastq_reads.values()) # total bp = sum of read lengths

    # write species_stats to output
    for k, v in species_sum.items():
        csv_templatestats.writerow([k] + v)

print(f"Unmapped fastq reads: {len(fastq_reads)}")

t3=time.perf_counter()
print(f"Processed {samlines} lines in {t3 -t2:0.4f} seconds")
