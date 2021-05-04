from __future__ import division
import sys
import argparse
import re
import os
import csv
import subprocess
from collections import OrderedDict
from collections import defaultdict
from collections import defaultdict
import glob
import readline
import pandas as pd
import timeit
import time
import gc
import datetime

# Parse Command line Arguments
parser = argparse.ArgumentParser(
    description='This script will parse StrainEst abundance results and generate an abundance matrix for the sample output in given input directory')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-StrainEst_dir', action='store', dest="StrainEst_dir",
                      help='StrainEst Input Directory - This directory contains StrainEst output results for the samples. Assumes the StrainEst output are placed in each sample name prefix directory. For ex: CDIF_SB_1206_/abund.txt')
args = parser.parse_args()

print "Getting StrainEst abundance output files - abund.txt in %s" % args.StrainEst_dir

abundance_files = glob.glob('%s/*/abund.txt' % args.StrainEst_dir)

# abundance_df = pd.DataFrame()


abundance_df = []
paste_command= "paste -d',' "

print "Reading first file %s to extract rownames" % abundance_files[0]
df = pd.read_csv("%s" % abundance_files[0], sep='\t', header=0)
#abundance_df.append(df[df.columns[0]])
df[df.columns[0]].to_csv('rownames', index=False, header=True)
paste_command= paste_command + "rownames "

for file in abundance_files:
    df1 = pd.read_csv("%s" % file, sep='\t', header=0)
    #abundance_df.append(df1[df1.columns[1]])
    filename = file.replace('abund.txt', 'abundance_values.txt')
    df1[df1.columns[1]].to_csv(filename, index=False, header=True)
    paste_command = paste_command + " %s" % filename

paste_command = paste_command + " > StrainEst_abundance_matrix.csv"
os.system(paste_command)

