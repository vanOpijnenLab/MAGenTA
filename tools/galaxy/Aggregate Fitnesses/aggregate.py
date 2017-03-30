# A translation of aggregate.pl into python! For analysis of Tn-Seq.
# This script requires BioPython just like calc_fitness.py, so you need it installed along with its dependencies if you want to run these scripts on your own.
# How to install BioPython and a list of its dependencies can be found here: http://biopython.org/DIST/docs/install/Installation.html
# K. McCoy









##### ARGUMENTS #####

def print_usage():
	print "Aggregate.py's usage is as follows:" + "\n\n"
	print "\033[1m" + "Required" + "\033[0m" + "\n"
	print "-o" + "\t\t" + "Output file for aggregated data." + "\n"
	print "\n"
	print "\033[1m" + "Optional" + "\033[0m" + "\n"
	print "-c" + "\t\t" + "Check for missing genes in the data set - provide a reference genome in genbank format. Missing genes will be sent to stdout." + "\n"
	print "-m" + "\t\t" + "Place a mark in an extra column for this set of genes. Provide a file with a list of genes seperated by newlines." + "\n"
	print "-x" + "\t\t" + "Cutoff: Don't include fitness scores with average counts (c1+c2)/2 < x (default: 0)" + "\n"
	print "-b" + "\t\t" + "Bottleneck value: The percentage of insertions randomly lost, which will be discounted for all genes (for example, 20% would be entered as 0.20; default 0.0)" + "\n"
	print "-f" + "\t\t" + "An in-between file carrying information on the blank count found from calc_fitness or consol_fitness; one of two ways to pass a blank count to this script" + "\n"
	print "-w" + "\t\t" + "Use weighted algorithm to calculate averages, variance, sd, se" + "\n"
	print "-l" + "\t\t" + "Weight ceiling: maximum value to use as a weight (default: 999,999)" + "\n"
	print "\n"
	print "All remainder arguements will be treated as fitness files (those files created by calc_fitness.py)" + "\n"
	print "\n"
	
import argparse 
parser = argparse.ArgumentParser()
parser.add_argument("-o", action="store", dest="summary")
parser.add_argument("-c", action="store", dest="find_missing")
parser.add_argument("-m", action="store", dest="marked")
parser.add_argument("-x", action="store", dest="cutoff")
parser.add_argument("-b", action="store", dest="blank_pc")
parser.add_argument("-f", action="store", dest="blank_file")
parser.add_argument("-w", action="store", dest="weighted")
parser.add_argument("-l", action="store", dest="weight_ceiling")
parser.add_argument("fitnessfiles", nargs=argparse.REMAINDER)

arguments = parser.parse_args()

if not arguments.summary:
	print "\n" + "You are missing a value for the -o flag. "
	print_usage() 
	quit()

if not arguments.fitnessfiles:
	print "\n" + "You are missing fitness file(s); these should be entered immediately after all the flags. "
	print_usage() 
	quit()
	
# 999,999 is a trivial placeholder number
	
if (not arguments.weight_ceiling):
	arguments.weight_ceiling = 999999
	
# Cutoff exists to discard positions with a low number of counted transcripts, because their fitness may not be as accurate - for the same reasoning that studies with low sample sizes can be innacurate. 
	
if (not arguments.cutoff):
	arguments.cutoff = 0

# Gets information from the txt output file of calc_fit / consol, if inputted

if arguments.blank_file:
	with open(arguments.blank_file) as file:
		blank_pc = file.read().splitlines()
		arguments.blank_pc = float(blank_pc[0].split()[1])

if (not arguments.blank_pc):
	arguments.blank_pc = 0





##### SUBROUTINES #####

# A subroutine that calculates the average, variance, standard deviation (sd), and standard error (se) of a group of scores; for use when aggregating scores by gene later on

import math
def unweighted_average(scores):
	sum = 0
	num = 0
	i = 0
	while i < len(scores):
		if not scores[i]:
			scores[i] = 0.0
		sum += float(scores[i])
		num += 1
		i += 1
	average = sum/num
	xminusxbars = 0
	while i < len(scores):
		xminusxbars += (float(scores[i]) - average)**2
	if num <= 1:
		variance = 0
	else:
		variance = xminusxbars/(num-1)
	sd = math.sqrt(variance)
	se = sd / math.sqrt(num)
	return (average, variance, sd, se)

# A subroutine that calculates the weighted average, variance, standard deviation (sd), and standard error (se) of a group of scores; the weights come from the number of reads each insertion location has
# For use when aggregating scores by gene later on, if the weighted argument is called

def weighted_average(scores,weights):
	sum = 0
	weighted_average = 0
	weighted_variance = 0
	top = 0
	bottom = 0
	i = 0
	while i < len(weights):
		if not scores[i]:
			scores[i] = 0.0
		top += float(weights[i])*float(scores[i])
		bottom += float(weights[i])
		i += 1
	if bottom == 0:
		return 0
	weighted_average = top/bottom
	top = 0
	bottom = 0
	i = 0
	while i < len(weights):
		top += float(weights[i]) * (float(scores[i]) - weighted_average)**2
		bottom += float(weights[i])
		i += 1
	weighted_variance = top/bottom
 	weighted_stdev = math.sqrt(weighted_variance)
	weighted_stder = weighted_stdev/math.sqrt(len(scores))
	return (weighted_average, weighted_variance, weighted_stdev, weighted_stder)










##### AGGREGATION / CALCULATIONS #####

#Reads the genes which should be marked in the final aggregate file into an array

import os.path
if arguments.marked:
	with open(arguments.marked) as file:
		marked_set = file.read().splitlines()

#Creates a dictionary of dictionaries to contain a summary of all genes and their fitness values
#The fitness values and weights match up, so that the weight of gene_summary[locus]["w"][2] would be gene_summary[locus]["s"][2]

import csv
gene_summary = {}
for eachfile in arguments.fitnessfiles:
	with open(eachfile) as csvfile:
		lines = csv.reader(csvfile)
		for line in lines:
			locus = line[9]
			w = line[12]
			if w == 'nW':
				continue
			if not w:
				w == 0
			c1 = float(line[2])
			c2 = float(line[3])
			avg = (c1+c2)/2
			if avg < float(arguments.cutoff):
				continue
			if avg > float(arguments.weight_ceiling):
				avg = arguments.weight_ceiling
			if locus not in gene_summary:
				gene_summary[locus] = {"w" : [], "s": []}
			gene_summary[locus]["w"].append(w)
			gene_summary[locus]["s"].append(avg)

#If finding any missing gene loci is requested in the arguments, starts out by loading all the known features from a genbank file

from Bio import SeqIO
if (arguments.find_missing):
	output = [["locus","mean","var","sd","se","gene","Total","Blank","Not Blank","Blank Removed","M\n"]]
	handle = open(arguments.find_missing, "rU")
	for record in SeqIO.parse(handle, "genbank"):
		refname = record.id
		features = record.features
	handle.close()
	
#Goes through the features to find which are genes
	
	for feature in features:
		gene = ""
		if feature.type == "gene":
			locus = "".join(feature.qualifiers["locus_tag"])
			if "gene" in feature.qualifiers:
				gene = "".join(feature.qualifiers["gene"])
		else:
			continue
			
#Goes through the fitness scores of insertions within each gene, and removes whatever % of blank fitness scores were requested along with their corresponding weights

		sum = 0
		num = 0
		avgsum = 0
		blank_ws = 0
		i = 0
		if locus in gene_summary.keys():
			for w in gene_summary[locus]["w"]:
				if float(w) == 0:
					blank_ws += 1
				else:
					sum += float(w)
					num += 1
			count = num + blank_ws		
			removed = 0
			to_remove = int(float(arguments.blank_pc)*count)
			if blank_ws > 0:
				i = 0
				while i < len(gene_summary[locus]["w"]):
					w = gene_summary[locus]["w"][i] 
					if removed == to_remove:
						break
					if float(w) == 0:
						del gene_summary[locus]["w"][i]
						del gene_summary[locus]["s"][i]
						removed += 1
						i -= 1
					i += 1
			
#If all the fitness values within a gene are empty, sets mean/var to 0.10 and Xs out sd/se; marks the gene if that's requested

			if num == 0:	
				if (arguments.marked and locus in marked_set):
					output.append([locus, "0.10", "0.10", "X", "X", gene, count, blank_ws, num, removed, "M", "\n"])
				else:
					output.append([locus, "0.10", "0.10", "X", "X", gene, count, blank_ws, num, removed, "\n"])

#Otherwise calls average() or weighted_average() to find the aggregate w / count / standard deviation / standard error of the insertions within each gene; marks the gene if that's requested

			else:
				if not arguments.weighted:
					(average, variance, stdev, stderr) = unweighted_average(gene_summary[locus]["w"])
				else:
					(average, variance, stdev, stderr) = weighted_average(gene_summary[locus]["w"],gene_summary[locus]["s"])
				if (arguments.marked and locus in marked_set):
					output.append([locus, average, variance, stdev, stderr, gene, count, blank_ws, num, removed, "M", "\n"])
				else:
					output.append([locus, average, variance, stdev, stderr, gene, count, blank_ws, num, removed, "\n"])
		
#If a gene doesn't have any insertions, sets mean/var to 0.10 and Xs out sd/se, plus leaves count through removed blank because there were no reads.

		else:
			if (arguments.marked and locus in marked_set):
				output.append([locus, "0.10", "0.10", "X", "X", gene, "", "", "", "", "M", "\n"])
			else:
				output.append([locus, "0.10", "0.10", "X", "X", gene, "", "", "", "", "\n"])

#Writes the aggregated fitness file

	with open(arguments.summary, "wb") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerows(output)

#If finding missing genes is not requested, just finds the aggregate w / count / standard deviation / standard error of the insertions within each gene, and writes them to a file, plus marks the genes requested
#This is never called through Galaxy since finding missing genes is just better than not finding them.

else:
	output = [["Locus","W","Count","SD","SE","M\n"]]
	for gene in gene_summary.keys():
		sum = 0
		num = 0
		average = 0
		if "w" not in gene_summary[gene]:
			continue
		for i in gene_summary[gene]["w"]:
			sum += i
			num += 1
		average = sum/num
		xminusxbars = 0
		for i in w:
			xminusxbars += (i-average)**2
		if num > 1:
			sd = math.sqrt(xminusxbars/(num-1))
			se = sd / math.sqrt(num)
		if (arguments.marked and locus in marked_set):
			output.append([gene, average, num, sd, se, "M", "\n"])
		else:
			output.append([gene, average, num, sd, se, "\n"])
	with open(arguments.summary, "wb") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerows(output)
