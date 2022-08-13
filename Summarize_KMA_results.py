import pandas as pd
import logging
import argparse
import re
import sys

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')
parser = argparse.ArgumentParser(description='Summarize KMA result files for one or more experiments')
parser.add_argument("input", nargs='+', metavar="--input", help="Filepath(s) of KMA result file(s) (.res with corresponding .mapstat)")
parser.add_argument("-t", metavar="--taxonomic_level", type=str, help="Taxonomic level at which to aggregate data, default: template", default="template")
parser.add_argument("-m", metavar="--experiment_metadata", type=str, help="Filepath to optional experiment metadata, excel ;-separated csv", default=None)
args = parser.parse_args()

# dictionary to couple .res filepath to .mapstat filepath
filepaths = {fp : re.sub(".res", ".mapstat", fp) for fp in args.input}

# read .res and .mapstat files, merge, and add Experiment variable
data = []
for res, mapstat in filepaths.items():
    experiment=re.search("\w+.res", res).group()
    KMAres = pd.read_csv(res, sep='\t')
    KMAmapstat = pd.read_csv(mapstat, sep='\t', skiprows=6)
    KMA = pd.merge(KMAmapstat, KMAres, how='outer', left_on='# refSequence', right_on='#Template')
    KMA["Experiment"] = str(experiment) # add column with experiment from which the data originates

    data.append(KMA)

# concatenate all data
KMAresults = pd.concat(data)

# optionally add metadata
if args.experiment_metadata!=None :
    metadata = pd.read_csv(args.experiment_metadata, sep=";")
    KMAresults = pd.merge(KMAresults, metadata, how='outer', on='Experiment')

# print summarized results to stdout in tsv format
KMAresults.to_csv(sys.stdout, sep="\t")