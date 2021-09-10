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
from collections import Counter
from argparse import RawTextHelpFormatter

# Parse Command line Arguments
parser = argparse.ArgumentParser(
    description='Generate Mixed Strain Simulation data for Strainest analysis.\n\n'
                'Using these four samples for various combinations: CDIF 840 (ST 110), CDIF 70 (ST 2), CDIF 399 (ST 8), CDIF 1000 (ST 231). \n\nUsing these two pairs for different relative abundance mixtures (99/1, 90/10, 80/20, 70/30, 60/40 and 50/50)\n\nCDIF 70 + CDIF 399\nCDIF 1000 + CDIF 840', formatter_class=RawTextHelpFormatter)
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-sample1', action='store', dest="sample1",
                      help='Sample 1 in the mixture')
required.add_argument('-sample2', action='store', dest="sample2",
                      help='Sample 2 in the mixture')
args = parser.parse_args()

# Get Total number of reads in Sample 1
proc = subprocess.Popen(["echo $(zcat %s|wc -l)/4|bc" % args.sample1], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()
total_reads_sample1 = int(out.strip())

# Get Total number of reads in Sample 2
proc = subprocess.Popen(["echo $(zcat %s|wc -l)/4|bc" % args.sample2], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()
total_reads_sample2 = int(out.strip())

print "\nTotal reads in Sample 1 - %s" % total_reads_sample1
print "\nTotal reads in Sample 2 - %s" % total_reads_sample2

# Iteration 1 where the base total is the total number of reads in Sample 1
print "\nIteration 1 where the base total is the total number of reads in Sample 1 - %s i.e %s" % (args.sample1, total_reads_sample1)
base_total = total_reads_sample1
base_sample = (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', '')


perc_99_basetotal = int(0.99 * base_total)
perc_1_basetotal = int(0.01 * base_total)
perc_90_basetotal = int(0.9 * base_total)
perc_10_basetotal = int(0.1 * base_total)
perc_80_basetotal = int(0.8 * base_total)
perc_20_basetotal = int(0.2 * base_total)
perc_70_basetotal = int(0.7 * base_total)
perc_30_basetotal = int(0.3 * base_total)
perc_60_basetotal = int(0.6 * base_total)
perc_40_basetotal = int(0.4 * base_total)
perc_50_basetotal = int(0.5 * base_total)


# Check the calculations
if perc_99_basetotal < base_total and perc_1_basetotal < total_reads_sample2:
    samplename_99_1 = "%s_99perc_%s_1perc" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 99/1: %s, %s, %s into %s" % (perc_99_basetotal, perc_1_basetotal, base_total, samplename_99_1)

    seqtk_base_sample = "seqtk sample -s100 %s %s > /tmp/%s_99perc_R1.fastq" % (args.sample1, perc_99_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s100 %s %s > /tmp/%s_1perc_R1.fastq" % (args.sample2, perc_1_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_99perc_R1.fastq /tmp/%s_1perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_99_1, samplename_99_1)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s100 %s %s > /tmp/%s_99perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_99_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s100 %s %s > /tmp/%s_1perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_1_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_99perc_R2.fastq /tmp/%s_1perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_99_1, samplename_99_1)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"

else:
    print "99/1 Not enough reads in Sample 2: %s < %s" % (perc_1_basetotal, total_reads_sample2)


if perc_90_basetotal < base_total and perc_10_basetotal < total_reads_sample2:
    samplename_90_10 = "%s_90perc_%s_10perc" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 90/10: %s, %s, %s into %s" % (perc_90_basetotal, perc_10_basetotal, base_total, samplename_90_10)

    seqtk_base_sample = "seqtk sample -s100 %s %s > /tmp/%s_90perc_R1.fastq" % (args.sample1, perc_90_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s100 %s %s > /tmp/%s_10perc_R1.fastq" % (args.sample2, perc_10_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_90perc_R1.fastq /tmp/%s_10perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_90_10, samplename_90_10)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s100 %s %s > /tmp/%s_90perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_90_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s100 %s %s > /tmp/%s_10perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_10_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_90perc_R2.fastq /tmp/%s_10perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_90_10, samplename_90_10)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"
else:
    print "90/10 Not enough reads in Sample 2: %s < %s" % (perc_10_basetotal, total_reads_sample2)

if perc_80_basetotal < base_total and perc_20_basetotal < total_reads_sample2:
    samplename_80_20 = "%s_80perc_%s_20perc" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 80/20: %s, %s, %s into %s" % (perc_80_basetotal, perc_20_basetotal, base_total, samplename_80_20)

    seqtk_base_sample = "seqtk sample -s200 %s %s > /tmp/%s_80perc_R1.fastq" % (args.sample1, perc_80_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s200 %s %s > /tmp/%s_20perc_R1.fastq" % (args.sample2, perc_20_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_80perc_R1.fastq /tmp/%s_20perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_80_20, samplename_80_20)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s200 %s %s > /tmp/%s_80perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_80_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s200 %s %s > /tmp/%s_20perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_20_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_80perc_R2.fastq /tmp/%s_20perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_80_20, samplename_80_20)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"
else:
    print "80/20 Not enough reads in Sample 2: %s < %s" % (perc_20_basetotal, total_reads_sample2)

if perc_70_basetotal < base_total and perc_30_basetotal < total_reads_sample2:
    samplename_70_30 = "%s_70perc_%s_30perc" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 70/30: %s, %s, %s into %s" % (perc_70_basetotal, perc_30_basetotal, base_total, samplename_70_30)

    seqtk_base_sample = "seqtk sample -s300 %s %s > /tmp/%s_70perc_R1.fastq" % (args.sample1, perc_70_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s300 %s %s > /tmp/%s_30perc_R1.fastq" % (args.sample2, perc_30_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_70perc_R1.fastq /tmp/%s_30perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_70_30, samplename_70_30)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s300 %s %s > /tmp/%s_70perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_70_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s300 %s %s > /tmp/%s_30perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_30_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_70perc_R2.fastq /tmp/%s_30perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_70_30, samplename_70_30)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"
else:
    print "70/30 Not enough reads in Sample 2: %s < %s" % (perc_30_basetotal, total_reads_sample2)

if perc_60_basetotal < base_total and perc_40_basetotal < total_reads_sample2:
    samplename_60_40 = "%s_60perc_%s_40perc" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 60/40: %s, %s, %s into %s" % (perc_60_basetotal, perc_40_basetotal, base_total, samplename_60_40)

    seqtk_base_sample = "seqtk sample -s400 %s %s > /tmp/%s_60perc_R1.fastq" % (args.sample1, perc_60_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s400 %s %s > /tmp/%s_40perc_R1.fastq" % (args.sample2, perc_40_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_60perc_R1.fastq /tmp/%s_40perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_60_40, samplename_60_40)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s400 %s %s > /tmp/%s_60perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_60_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s400 %s %s > /tmp/%s_40perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_40_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_60perc_R2.fastq /tmp/%s_40perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_60_40, samplename_60_40)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"
else:
    print "60/40 Not enough reads in Sample 2: %s < %s" % (perc_40_basetotal, total_reads_sample2)

if perc_50_basetotal < base_total and perc_50_basetotal < total_reads_sample2:
    samplename_50_50 = "%s_50perc_%s_50perc" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 50/50: %s, %s, %s into %s" % (perc_50_basetotal, perc_50_basetotal, base_total, samplename_50_50)

    seqtk_base_sample = "seqtk sample -s500 %s %s > /tmp/%s_50perc_R1.fastq" % (args.sample1, perc_50_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s500 %s %s > /tmp/%s_50perc_R1.fastq" % (args.sample2, perc_50_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_50perc_R1.fastq /tmp/%s_50perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_50_50, samplename_50_50)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s500 %s %s > /tmp/%s_50perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_50_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s500 %s %s > /tmp/%s_50perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_50_basetotal, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_50perc_R2.fastq /tmp/%s_50perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', ''), samplename_50_50, samplename_50_50)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"

else:
    print "50/50 Not enough reads in Sample 2: %s < %s" % (perc_50_basetotal, total_reads_sample2)

# Iteration 1 where the base total is the total number of reads in Sample 1
print "\nIteration 2 where the base total is the total number of reads in Sample 2 - %s i.e %s" % (args.sample2, total_reads_sample2)
base_total = total_reads_sample2
base_sample = (os.path.basename(args.sample2)).replace('_R1_001.fastq.gz', '')


perc_99_basetotal = int(0.99 * base_total)
perc_1_basetotal = int(0.01 * base_total)
perc_90_basetotal = int(0.9 * base_total)
perc_10_basetotal = int(0.1 * base_total)
perc_80_basetotal = int(0.8 * base_total)
perc_20_basetotal = int(0.2 * base_total)
perc_70_basetotal = int(0.7 * base_total)
perc_30_basetotal = int(0.3 * base_total)
perc_60_basetotal = int(0.6 * base_total)
perc_40_basetotal = int(0.4 * base_total)
perc_50_basetotal = int(0.5 * base_total)


# Check the calculations
if perc_99_basetotal < base_total and perc_1_basetotal < total_reads_sample1:
    samplename_99_1 = "%s_99perc_%s_1perc" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 99/1: %s, %s, %s into %s" % (perc_99_basetotal, perc_1_basetotal, base_total, samplename_99_1)

    seqtk_base_sample = "seqtk sample -s100 %s %s > /tmp/%s_99perc_R1.fastq" % (args.sample2, perc_99_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s100 %s %s > /tmp/%s_1perc_R1.fastq" % (args.sample1, perc_1_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_99perc_R1.fastq /tmp/%s_1perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_99_1, samplename_99_1)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s100 %s %s > /tmp/%s_99perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_99_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s100 %s %s > /tmp/%s_1perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_1_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_99perc_R2.fastq /tmp/%s_1perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_99_1, samplename_99_1)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"

else:
    print "99/1 Not enough reads in Sample 1: %s < %s" % (perc_1_basetotal, total_reads_sample1)

if perc_90_basetotal < base_total and perc_10_basetotal < total_reads_sample1:
    samplename_90_10 = "%s_90perc_%s_10perc" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 90/10: %s, %s, %s into %s" % (perc_90_basetotal, perc_10_basetotal, base_total, samplename_90_10)

    seqtk_base_sample = "seqtk sample -s100 %s %s > /tmp/%s_90perc_R1.fastq" % (args.sample2, perc_90_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s100 %s %s > /tmp/%s_10perc_R1.fastq" % (args.sample1, perc_10_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_90perc_R1.fastq /tmp/%s_10perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_90_10, samplename_90_10)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s100 %s %s > /tmp/%s_90perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_90_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s100 %s %s > /tmp/%s_10perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_10_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_90perc_R2.fastq /tmp/%s_10perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_90_10, samplename_90_10)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"

else:
    print "90/10 Not enough reads in Sample 1: %s < %s" % (perc_10_basetotal, total_reads_sample1)

if perc_80_basetotal < base_total and perc_20_basetotal < total_reads_sample1:
    samplename_80_20 = "%s_80perc_%s_20perc" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 80/20: %s, %s, %s into %s" % (perc_80_basetotal, perc_20_basetotal, base_total, samplename_80_20)

    seqtk_base_sample = "seqtk sample -s200 %s %s > /tmp/%s_80perc_R1.fastq" % (args.sample2, perc_80_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s200 %s %s > /tmp/%s_20perc_R1.fastq" % (args.sample1, perc_20_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_80perc_R1.fastq /tmp/%s_20perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_80_20, samplename_80_20)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s200 %s %s > /tmp/%s_80perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_80_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s200 %s %s > /tmp/%s_20perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_20_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_80perc_R2.fastq /tmp/%s_20perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_80_20, samplename_80_20)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"
else:
    print "80/20 Not enough reads in Sample 1: %s < %s" % (perc_20_basetotal, total_reads_sample1)

if perc_70_basetotal < base_total and perc_30_basetotal < total_reads_sample1:
    samplename_70_30 = "%s_70perc_%s_30perc" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 70/30: %s, %s, %s into %s" % (perc_70_basetotal, perc_30_basetotal, base_total, samplename_70_30)

    seqtk_base_sample = "seqtk sample -s300 %s %s > /tmp/%s_70perc_R1.fastq" % (args.sample2, perc_70_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s300 %s %s > /tmp/%s_30perc_R1.fastq" % (args.sample1, perc_30_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_70perc_R1.fastq /tmp/%s_30perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_70_30, samplename_70_30)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s300 %s %s > /tmp/%s_70perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_70_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s300 %s %s > /tmp/%s_30perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_30_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_70perc_R2.fastq /tmp/%s_30perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_70_30, samplename_70_30)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"
else:
    print "70/30 Not enough reads in Sample 1: %s < %s" % (perc_30_basetotal, total_reads_sample1)

if perc_60_basetotal < base_total and perc_40_basetotal < total_reads_sample1:
    samplename_60_40 = "%s_60perc_%s_40perc" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 60/40: %s, %s, %s into %s" % (perc_60_basetotal, perc_40_basetotal, base_total, samplename_60_40)

    seqtk_base_sample = "seqtk sample -s400 %s %s > /tmp/%s_60perc_R1.fastq" % (args.sample2, perc_60_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s400 %s %s > /tmp/%s_40perc_R1.fastq" % (args.sample1, perc_40_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_60perc_R1.fastq /tmp/%s_40perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_60_40, samplename_60_40)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s400 %s %s > /tmp/%s_60perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_60_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s400 %s %s > /tmp/%s_40perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_40_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_60perc_R2.fastq /tmp/%s_40perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_60_40, samplename_60_40)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"

else:
    print "60/40 Not enough reads in Sample 1: %s < %s" % (perc_40_basetotal, total_reads_sample1)

if perc_50_basetotal < base_total and perc_50_basetotal < total_reads_sample1:
    samplename_50_50 = "%s_50perc_%s_50perc" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    print "#####################\nSimulating 50/50: %s, %s, %s into %s" % (perc_50_basetotal, perc_50_basetotal, base_total, samplename_50_50)

    seqtk_base_sample = "seqtk sample -s500 %s %s > /tmp/%s_50perc_R1.fastq" % (args.sample2, perc_50_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s500 %s %s > /tmp/%s_50perc_R1.fastq" % (args.sample1, perc_50_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_50perc_R1.fastq /tmp/%s_50perc_R1.fastq > %s_R1.fastq && gzip %s_R1.fastq" % (base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_50_50, samplename_50_50)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)

    seqtk_base_sample = "seqtk sample -s500 %s %s > /tmp/%s_50perc_R2.fastq" % (
        (args.sample2).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_50_basetotal, base_sample)
    seqtk_sample2 = "seqtk sample -s500 %s %s > /tmp/%s_50perc_R2.fastq" % (
        (args.sample1).replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'), perc_50_basetotal, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''))
    combine_samples = "cat /tmp/%s_50perc_R2.fastq /tmp/%s_50perc_R2.fastq > %s_R2.fastq && gzip %s_R2.fastq" % (
    base_sample, (os.path.basename(args.sample1)).replace('_R1_001.fastq.gz', ''), samplename_50_50, samplename_50_50)

    os.system(seqtk_base_sample)
    os.system(seqtk_sample2)
    os.system(combine_samples)
    print "#####################\n"

else:
    print "50/50 Not enough reads in Sample 1: %s < %s" % (perc_50_basetotal, total_reads_sample1)

