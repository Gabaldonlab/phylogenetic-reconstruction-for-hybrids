#!/usr/bin/env python

#Build strain tree

import sys
sys.path.append("/path/to/your/scripts/")
import genome_modules_3 as GM
from tqdm import tqdm
import glob
import argparse

#Detect positions that should be omitted from the alignment because they are badly covered. It assumes that you have:

# 1.- a reference genome
# 2.- a folder with all the bed files resulting from the genomecov analysis
# 3.- a folder with the filtered vcf files.
# 4.- Strain names will be taken from fileNames in folder 3


parser = argparse.ArgumentParser(description="Create pseudo alignment")
parser.add_argument("-r",dest="referenceFile",action="store",required=True,help="Reference genome")
parser.add_argument("-c",dest="coverageFolder",action="store",required=True,help="Folder where genome cov results are saved")
parser.add_argument("-s",dest="snpFolder",action="store",required=True,help="Folder where filtered SNPs are saved")
parser.add_argument("-o",dest="outFile",action="store",required=True,help="Outfile name")
parser.add_argument("--min_reads",dest="min_reads",action="store",type=int,default=5,help="Minimum number of reads so that a position is kept")
args = parser.parse_args()

#Load the reference genome. This will serve as base to build the alignment
try:
    seqs = GM.load_sequences(args.referenceFile," ")
except:
    exit("unable to load the reference genome")

positions = {}
#Scan genomecov file to detect positions that are badly covered
for genomeCov_file in tqdm(glob.glob(args.coverageFolder+"/*")):
    total_coverage = {}
    for line in open(genomeCov_file):
        line = line.strip()
        dades = line.split("\t")
        if int(dades[-1]) < args.min_reads:
            if dades[0] not in positions:
                positions[dades[0]] = set([])
            if dades[0] not in total_coverage:
                total_coverage[dades[0]] = 0
            for a in range(int(dades[1]),int(dades[2])):
                positions[dades[0]].add(a)
                total_coverage[dades[0]] += 1
    # ~ for chromo in total_coverage:
        # ~ p = total_coverage[chromo] / len(seqs[chromo]) *100.0
        # ~ print(strain,chromo,total_coverage[chromo],len(seqs[chromo]),p)

#Load SNP data
snps = {}
positions_with_snps = {}
strains = set([])
for fileName in tqdm(glob.glob(args.snpFolder+"/*")):
    strain = fileName.split("/")[-1].split(".")[0]
    strains.add(strain)
    snps[strain] = {}
    for line in open(fileName):
        line = line.strip()
        dades = line.split("\t")
        if "#" in line:
            pass
        elif "PASS" in line:
            dades = line.split("\t")
            pos = int(dades[1]) - 1
            ref = dades[3]
            alt = dades[4]
            contig = dades[0]
            if pos not in positions[contig]:
                #If there's an indel or a deletion, omit the position as you can't be sure what to do with it
                if len(ref) > 1 or len(alt) > 1:
                    positions[contig].add(pos)
                else:
                    if contig not in snps[strain]:
                        snps[strain][contig] = {}
                    if contig not in positions_with_snps:
                        positions_with_snps[contig] = set([])
                    positions_with_snps[contig].add(pos)
                    if "1/1" in line:
                        snps[strain][contig][pos] = [ref,alt] # for homozygous SNPs take the alternative allele
                    elif "0/1" in line:
                        snps[strain][contig][pos] = [ref,alt] # for heterozygous SNPs take the alternative allele only

#Build a file for each chromosome
for contig in tqdm(seqs):
    outfile = open(contig+".fa","w")
    reference = list(seqs[contig])
    for strain in strains:
        string = ""
        changes_made = 0
        for a in range(len(seqs[contig])):
            if contig in positions and contig in positions_with_snps:
                if a not in positions[contig]:
                    if a in positions_with_snps[contig]:
                        if contig not in snps[strain]:
                            string += reference[a]
                        else:
                            if a not in snps[strain][contig]:
                                string += reference[a]
                            else:
                                if snps[strain][contig][a][0] != reference[a]:
                                    print("CHECK POS",strain,contig,a,snps[strain][contig][a])
                                string += snps[strain][contig][a][1] # here is where the script adds the snp, [1] refers to position 1 in the snps dicctionary in lines 83 and 85
                                changes_made += 1
        print("CHANGES:",strain,contig,changes_made)
        GM.print_sequence(strain,string,outfile)
    string = ""
    for a in range(len(seqs[contig])):
        if contig in positions and contig in positions_with_snps:
            if a not in positions[contig]:
                if a in positions_with_snps[contig]:
                    string += reference[a]
    GM.print_sequence("Reference",string,outfile)
    outfile.close()

#Concatenate alignments
concatenated = {}
for contig in seqs:
    seqs1 = GM.load_sequences(contig+".fa"," ")
    for strain in seqs1:
        if strain not in concatenated:
            concatenated[strain] = ""
        concatenated[strain] += seqs1[strain]

outfile = open(args.outFile,"w")
for strain in concatenated:
    GM.print_sequence(strain,concatenated[strain],outfile)
outfile.close()
