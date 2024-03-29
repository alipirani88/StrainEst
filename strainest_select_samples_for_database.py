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
import numpy as np
import timeit
import time
import gc
import datetime
from collections import Counter
from pyfasta import Fasta

# Parse Command line Arguments
parser = argparse.ArgumentParser(
    description='This script will parse Ariba MLST report, select a single ST representative sample from each ST group based on highest GATK %_bases_above_5 value, generate a strainest compatible SNP matrix for database generation')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-mlst_report', action='store', dest="mlst_report",
                      help='Ariba MLST Report - First column=samplename, second column=ST')
optional.add_argument('-depth_report', action='store', dest="depth_report",
                      help='GATK Depth of coverage report (generated by concatenating)')
optional.add_argument('-matrix', action='store', dest="matrix",
                      help='SNP allele Matrix - All samples')
optional.add_argument('-code_matrix', action='store', dest="code_matrix",
                      help='SNP code Matrix - All samples')
optional.add_argument('-reference', action='store', dest="reference",
                      help='Reference Genome Fasta file')
args = parser.parse_args()

def extract_mlst():
    print "Extract ST for each samples and generate a Sample -> ST map\n"
    mlst_map = dict()
    with open("%s" % args.mlst_report, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        next(csv_reader, None)
        for row in csv_reader:
            if row[1] in mlst_map:
                mlst_map[row[1]].append(row[0])
            else:
                mlst_map[row[1]] = [row[0]]
    return mlst_map

# mlst_map = extract_mlst()

def extract_depth():
    print "Extract percentage of bases_above_5 for each samples from %s and generating Sample to GATK Depth map\n" % args.depth_report
    depth_map = dict()
    with open("%s" % args.depth_report, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        next(csv_reader, None)
        for row in csv_reader:
            depth_map[row[0]] = row[7]
    return depth_map
# depth_map = extract_depth()

def select_ST_based_representative_sample():
    print "Extracting Representative Sample for each ST based on the Max Read Depth and writing it to - Selected_ST_representative_samples.txt\n"
    representative_sample = []
    fp = open("Selected_ST_representative_samples.txt", 'w+')
    for key in mlst_map:
        if len(mlst_map[key]) == 1:
            representative_sample.append(str(mlst_map[key][0]) + '_R1_001.fastq.gz')
            fp.write("ST%s\t%s\t%s\t%s" % (key, len(mlst_map[key]), str(mlst_map[key][0]), depth_map[str(mlst_map[key][0])]))
            #print "ST%s\t%s\t%s\t%s" % (key, len(mlst_map[key]), str(mlst_map[key][0]), depth_map[str(mlst_map[key][0])])
        else:
            depth = 0
            rep = ""
            for i in mlst_map[key]:
                new_depth = float(depth_map[i])
                #print new_depth
                if new_depth > depth:
                    depth = new_depth
                    rep = i
            #print "ST%s\t%s\t%s\t%s" % (key, len(mlst_map[key]), rep, depth_map[rep])
            fp.write("ST%s\t%s\t%s\t%s" % (key, len(mlst_map[key]), rep, depth_map[rep]))
            representative_sample.append(str(rep) + "_R1_001.fastq.gz")
    fp.close()
    return representative_sample

# representative_sample = select_ST_based_representative_sample()

# Print SNP allele Matrix of Representative sample to Selected_ST_representative_matrix.csv
def extract_representative_samples():
    print "Reading SNP matrix and subsetting %s to contain ST based Representative Samples - Selected_ST_representative_matrix.csv\n" % args.matrix
    global df
    df = pd.read_csv("%s" % args.matrix, sep=',', header=0)
    df1 = df[representative_sample]
    df1.to_csv('Selected_ST_representative_matrix.csv', index=False)
    del df1
    gc.collect()

# print "Time taken to read SNP matrix and extract Representative Samples: %s" % (timeit.timeit(extract_representative_samples, number=1))

def extract_matrix_positions():
    #Extract Position and Reference allele from SNP matrix
    print "Extracting Position from SNP matrix - %s to POS\n" % args.matrix
    global df
    df = pd.read_csv("%s" % args.matrix, sep=',', header=0)
    global pos
    pos = []
    for i in range(len(df)) :
      stringarray = (df.iloc[i, 0]).split(' ')
      pos.append(stringarray[3])
    #Print positions to POS file
    fileObj = open("POS", 'w+')
    fileObj.write('Pos' + '\n')
    for item in pos:
        fileObj.write(item + '\n')
    fileObj.close()
    return pos

# pos = extract_matrix_positions()

def extract_reference_allele():
    print "Extracting Reference Allele from Reference Fasta file - %s to REF\n" % args.reference
    # Get reference genome ID from reference fasta file
    get_reference = Fasta(args.reference)
    if len(get_reference.keys()) == 1:
        ref_id = get_reference.keys()
    print "The reference genome ID from reference genome - %s" % ref_id

    fileObj = open("REF", 'w+')
    fileObj.write('Ref' + '\n')
    for item in pos:
        ref_allele = str(get_reference.sequence({'chr': str(get_reference.keys()[0]), 'start': int(item), 'stop': int(item)}))
        fileObj.write(ref_allele + '\n')
    fileObj.close()

# extract_reference_allele()


# # Generate StrainEst Database matrix
# print "Generate Strainest Database matrix by pasting files(POS, REF, Selected_ST_representative_matrix.csv) to Strainest_Selected_ST_representative_db.csv\n"
# os.system("paste -d ',' POS REF Selected_ST_representative_matrix.csv > Strainest_Selected_ST_representative_db.csv")
#
# print "Adding Row annotation column to representative_sample array so that it can also be subset to Selected_ST_representative_code_matrix.csv\n"
# representative_sample.insert(0, "Type of SNP at POS > ALT functional=PHAGE_REPEAT_MASK locus_tag=locus_id strand=strand; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos|gene_symbol|product")


def extract_representative_samples_code_matrix():
    print "Reading SNP code matrix and extracting Representative Samples to Selected_ST_representative_code_matrix.csv.\n"
    global df
    df = pd.read_csv("%s" % args.code_matrix, sep=',', header=0)
    df1 = df[representative_sample]
    df1.to_csv('Selected_ST_representative_code_matrix.csv', index=False)
    del df1
    gc.collect()

# extract_representative_samples_code_matrix()
#print "Time taken to read SNP code matrix and extract Representative Samples: %s\n" % (timeit.timeit(extract_representative_samples_code_matrix, number=1))


def check_Selected_ST_representative_code_matrix_code(codematrix):
    print "Calculating distribution of SNP codes in Selected_ST_representative_code_matrix_core_db.csv\n"
    df = pd.read_csv("%s" % codematrix, sep='\t', header=0)
    fp = open("Selected_ST_representative_code_matrix_core_db_bargraph.csv", 'w+')
    fp.write("Sample, Reference Allele, Core Variant(1), Hard Filtered Variant(2), True Variant but filtered in another sample(3), SNP proximate to Indel(4), Unmapped(-1), Phage(-2), FQ only filtered(-3), MQ only filtered(-4)")
    for (columnName, columnData) in df.iteritems():
        n = len(columnData.values)
        core_variants = 0
        unmapped = 0
        ref_allele = 0
        phage = 0
        fq_pos = 0
        mq_pos = 0
        hard_filter = 0
        filtered_in_another_sample = 0
        indel_proximate = 0

        for i in range(n):
            if columnData[i] == 1:
                core_variants += 1

            if columnData[i] == -1:
                unmapped += 1

            if columnData[i] == 0:
                ref_allele += 1

            if columnData[i] == -2:
                phage += 1

            if columnData[i] == -3:
                fq_pos += 1

            if columnData[i] == -4:
                mq_pos += 1

            if columnData[i] == 2:
                hard_filter += 1

            if columnData[i] == 3:
                core_variants += 1

            if columnData[i] == 4:
                indel_proximate += 1
        fp.write("%s, %s, %s, %s, %s, %s, %s, %s, %s\n" % (columnName, ref_allele, core_variants, hard_filter, indel_proximate, unmapped, phage, fq_pos, mq_pos))
        #print "%s, %s, %s, %s, %s, %s, %s, %s, %s" % (columnName, ref_allele, core_variants, hard_filter, indel_proximate, unmapped, phage, fq_pos, mq_pos)
    fp.close()

def extract_core_positions_from_SNP_code_matrix():
    print "Reading SNP code matrix and extracting core variants.\n"
    position_label = OrderedDict()
    non_core_codes = ['2', '4', '-2', '-3', '-4', '-1']
    core_codes = ['0', '1', '3']
    variant_code = ['1', '3']
    core_positions = []
    with open("%s" % args.code_matrix, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')

        next(csv_reader, None)
        for row in csv_reader:
            position_label[row[0]] = row[1:]
    csv_file.close()
    count = 0
    for key in position_label:
        if set(core_codes) & set(position_label[key]):
            if set(non_core_codes) & set(position_label[key]):
                continue
            else:
                count += 1
                core_positions.append(key.split(' ')[3])

    print "Reading Strainest_Selected_ST_representative_db matrix and extracting core positions to Strainest_Selected_ST_representative_core_db.csv.\n"
    global df
    df = pd.read_csv("Strainest_Selected_ST_representative_db.csv", sep=',', header=0)
    df1 = df.loc[df['Pos'].isin(core_positions)]
    df1.to_csv('Strainest_Selected_ST_representative_core_db.csv', index=False)
    del df1
    gc.collect()

    fp = open("Selected_ST_representative_code_matrix_core_db.csv", 'w+')
    with open("Selected_ST_representative_code_matrix.csv", 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            # print row
            position = (row[0]).split(' ')[3]
            if position in core_positions or (row[0]).startswith('Type'):
                s = '\t'.join(row)
                fp.write(str(s) + '\n')
    fp.close()
    check_Selected_ST_representative_code_matrix_code("Selected_ST_representative_code_matrix_core_db.csv")

#extract_core_positions_from_SNP_code_matrix()


def extract_core_and_unmapped_positions_from_SNP_code_matrix():
    print "Reading SNP code matrix and extracting core/unmapped variants..."
    position_label = OrderedDict()
    non_core_codes = ['2', '4', '-2', '-3', '-4']
    core_unmapped_codes = ['0', '1', '-1', '3']
    core_unmapped_positions = []
    with open("Selected_ST_representative_code_matrix.csv", 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')

        next(csv_reader, None)
        for row in csv_reader:
            position_label[row[0]] = row[1:]
    csv_file.close()
    count = 0
    for key in position_label:
        if set(core_unmapped_codes) & set(position_label[key]):
            if set(non_core_codes) & set(position_label[key]):
                continue
            else:
                count += 1
                core_unmapped_positions.append(key.split(' ')[3])
    #print core_positions
    print "Reading Strainest_Selected_ST_representative_db matrix and extracting core/unmapped positions to Strainest_Selected_ST_representative_core_unmapped_db.csv..."
    global df
    df = pd.read_csv("Strainest_Selected_ST_representative_db.csv", sep=',', header=0)
    df1 = df.loc[df['Pos'].isin(core_unmapped_positions)]
    df1.to_csv('Strainest_Selected_ST_representative_core_unmapped_db.csv', index=False)
    del df1
    gc.collect()

    fp = open("Selected_ST_representative_code_matrix_core_unmapped_db.csv", 'w+')
    with open("Selected_ST_representative_code_matrix.csv", 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            # print row
            position = (row[0]).split(' ')[3]
            if position in core_unmapped_positions or (row[0]).startswith('Type'):
                s = '\t'.join(row)
                fp.write(str(s) + '\n')
    fp.close()
    check_Selected_ST_representative_code_matrix_code("Selected_ST_representative_code_matrix_core_unmapped_db.csv")

#extract_core_and_unmapped_positions_from_SNP_code_matrix()

def extract_core_and_unmapped_positions_from_SNP_code_matrix_for_MLST_genes():
    print "Extract Positions for MLST genes..."
    mlst_pos = []
    os.system("grep -Ew 'adk|atpA|dxr|glyA|recA|sodA|tpi' %s | cut -d' ' -f4 > MLST_gene_positions.txt" % args.code_matrix)
    with open("MLST_gene_positions.txt", 'r') as fpp:
        for lines in fpp:
            mlst_pos.append(str(lines.strip()))


    print "Reading SNP code matrix and extracting core/unmapped variants..."
    position_label = OrderedDict()
    non_core_codes = ['2', '4', '-2', '-3', '-4']
    core_unmapped_codes = ['0', '1', '-1', '3']
    core_unmapped_positions = []
    with open("Selected_ST_representative_code_matrix.csv", 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')

        next(csv_reader, None)
        for row in csv_reader:
            position_label[row[0]] = row[1:]
    csv_file.close()
    count = 0
    for key in position_label:
        if set(core_unmapped_codes) & set(position_label[key]):
            if set(non_core_codes) & set(position_label[key]):
                continue
            else:
                count += 1
                if str(key.split(' ')[3]) in mlst_pos:
                    core_unmapped_positions.append(key.split(' ')[3])
    print "Reading Strainest_Selected_ST_representative_db matrix and extracting core/unmapped positions to Strainest_Selected_ST_representative_core_unmapped_db.csv..."
    global df
    df = pd.read_csv("Strainest_Selected_ST_representative_db.csv", sep=',', header=0)
    df1 = df.loc[df['Pos'].isin(core_unmapped_positions)]
    df1.to_csv('Strainest_Selected_ST_representative_core_unmapped_db.csv', index=False)
    del df1
    gc.collect()

    fp = open("Selected_ST_representative_code_matrix_core_unmapped_db.csv", 'w+')
    with open("Selected_ST_representative_code_matrix.csv", 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            # print row
            position = (row[0]).split(' ')[3]
            if position in core_unmapped_positions or (row[0]).startswith('Type'):
                s = '\t'.join(row)
                fp.write(str(s) + '\n')
    fp.close()
    check_Selected_ST_representative_code_matrix_code("Selected_ST_representative_code_matrix_core_unmapped_db.csv")


# extract_core_and_unmapped_positions_from_SNP_code_matrix_for_MLST_genes()

#check_Selected_ST_representative_code_matrix_code()

#########################################################

def QC_passed_all_genomes_core():
    pos = extract_matrix_positions()
    extract_reference_allele()

    global df
    df = pd.read_csv("%s" % args.matrix, sep=',', header=0)
    global representative_sample
    representative_sample = list(df.columns[1:])

    df1 = df[representative_sample]
    df1.to_csv('Selected_ST_representative_matrix.csv', index=False)
    del df1
    gc.collect()

    print "Generate Strainest Database matrix by pasting files(POS, REF, Selected_ST_representative_matrix.csv) to Strainest_Selected_ST_representative_db.csv\n"
    os.system("paste -d ',' POS REF Selected_ST_representative_matrix.csv > Strainest_Selected_ST_representative_db.csv")

    print "Adding Row annotation column to representative_sample array so that it can also be subset to Selected_ST_representative_code_matrix.csv\n"
    representative_sample.insert(0, "Type of SNP at POS > ALT functional=PHAGE_REPEAT_MASK locus_tag=locus_id strand=strand; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos|gene_symbol|product")

    extract_representative_samples_code_matrix()

    extract_core_positions_from_SNP_code_matrix()


def QC_passed_all_genomes_core_and_unmapped():
    pos = extract_matrix_positions()
    extract_reference_allele()

    global df
    df = pd.read_csv("%s" % args.matrix, sep=',', header=0)
    global representative_sample
    representative_sample = list(df.columns[1:])

    df1 = df[representative_sample]
    df1.to_csv('Selected_ST_representative_matrix.csv', index=False)
    del df1
    gc.collect()

    print "Generate Strainest Database matrix by pasting files(POS, REF, Selected_ST_representative_matrix.csv) to Strainest_Selected_ST_representative_db.csv\n"
    os.system("paste -d ',' POS REF Selected_ST_representative_matrix.csv > Strainest_Selected_ST_representative_db.csv")

    print "Adding Row annotation column to representative_sample array so that it can also be subset to Selected_ST_representative_code_matrix.csv\n"
    representative_sample.insert(0, "Type of SNP at POS > ALT functional=PHAGE_REPEAT_MASK locus_tag=locus_id strand=strand; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos|gene_symbol|product")

    extract_representative_samples_code_matrix()

    extract_core_and_unmapped_positions_from_SNP_code_matrix()


def QC_passed_all_genomes_core_exclude_test_samples():
    test_samples = []
    with open("test_samples.txt") as fp:
        for line in fp:
            line = line.strip()
            test_samples.append(line)
    fp.close()


    global df
    df = pd.read_csv("%s" % args.matrix, sep=',', header=0)
    global representative_sample
    representative_sample = list(df.columns[0:])
    print len(representative_sample)
    for i in representative_sample:
        if i in test_samples:
            representative_sample.remove(i)
    print len(representative_sample)
    df1 = df[representative_sample]
    df1.to_csv('Strainest_Selected_ST_representative_core_db_exclude_test_samples.csv', index=False)
    del df1
    gc.collect()

def QC_passed_all_genomes_core_and_unmapped_exclude_test_samples():
    test_samples = []
    with open("test_samples.txt") as fp:
        for line in fp:
            line = line.strip()
            test_samples.append(line)
    fp.close()


    global df
    df = pd.read_csv("%s" % args.matrix, sep=',', header=0)
    global representative_sample
    representative_sample = list(df.columns[0:])
    print len(representative_sample)
    for i in representative_sample:
        if i in test_samples:
            representative_sample.remove(i)
    print len(representative_sample)
    df1 = df[representative_sample]
    df1.to_csv('Strainest_Selected_ST_representative_core_db_exclude_test_samples.csv', index=False)
    del df1
    gc.collect()

def subset_strainest_db():
    print "Reading SNP code matrix and extracting core variants.\n"

    fp = open("Strainest_Selected_ST_representative_core_db_exclude_only_reference.csv", 'w+')
    with open("%s" % args.matrix, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            fp.write(','.join(row) + '\n')
            break

        next(csv_reader, None)
        for row in csv_reader:
            result = False
            if len(row) > 0:
                result = row.count(row[1]) == len(row[1:])
            if result:
                print row[0]
            else:
                fp.write(','.join(row) + '\n')
    csv_file.close()

    fp.close()

#########################################################
#mlst_map = extract_mlst()

#depth_map = extract_depth()

#representative_sample = select_ST_based_representative_sample()

#print "Time taken to read SNP matrix and extract Representative Samples: %s" % (timeit.timeit(extract_representative_samples, number=1))

# pos = extract_matrix_positions()

#extract_reference_allele()

# Generate StrainEst Database matrix
# print "Generate Strainest Database matrix by pasting files(POS, REF, Selected_ST_representative_matrix.csv) to Strainest_Selected_ST_representative_db.csv\n"
# os.system("paste -d ',' POS REF Selected_ST_representative_matrix.csv > Strainest_Selected_ST_representative_db.csv")
#
# print "Adding Row annotation column to representative_sample array so that it can also be subset to Selected_ST_representative_code_matrix.csv\n"
# representative_sample.insert(0, "Type of SNP at POS > ALT functional=PHAGE_REPEAT_MASK locus_tag=locus_id strand=strand; ALT|Effect|Impact|GeneID|Nrchange|Aachange|Nrgenepos|AAgenepos|gene_symbol|product")
#
# extract_representative_samples_code_matrix()
#
# extract_core_positions_from_SNP_code_matrix()
#
# extract_core_and_unmapped_positions_from_SNP_code_matrix()
#
# extract_core_and_unmapped_positions_from_SNP_code_matrix_for_MLST_genes()

#QC_passed_all_genomes_core()

#QC_passed_all_genomes_core_and_unmapped()

#QC_passed_all_genomes_core_exclude_test_samples()

QC_passed_all_genomes_core_and_unmapped_exclude_test_samples()

#subset_strainest_db()