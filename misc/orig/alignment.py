#A file to align fitness / aggregate type files by gene or insertion loci
#K McCoy yo

import argparse
import csv

def print_usage():
	print "Alignment.py's usage is as follows (all flags required):" + "\n\n"
	print "-r" + "\t\t" + "the reference file to which all other files will be aligned" + "\n"
	print "-o" + "\t\t" + "the name of the output file" + "\n"
	print "-c" + "\t\t" + "the columns of the reference file to be included in the new file, counting from 0, separated by commas" + "\n"
	print "-ac" + "\t\t" + "the column of the reference file used for alignment, counting from 0" + "\n"
	print "-u" + "\t\t" + "print the unaligned files to the output file, below the aligned rows" + "\n"
	print "-con" + "\t\t" + "a conversion file, if you're aligning different strains. Each line should have a gene from the refile with all corresponding genes following separated by tabs, like so: SPA0001	SPT0001	SPD0001" + "\n"
	print "\n"
	print "All remainder arguements should be separated by colons with the file name as the first value, the columns to be included in the new file as the second, and the column to align by as the third, counting from 0, like so: example.txt:0,3,4:0. These files should be excel/csv type files, and have header values in the first row and values to be aligned in the rest." + "\n"
	print "\n"


parser = argparse.ArgumentParser()
parser.add_argument("-r", action="store", dest="referencefile")
parser.add_argument("-o", action="store", dest="output")
parser.add_argument("-c", action="store", dest="refcolumns")
parser.add_argument("-ac", action="store", dest="refalignmentcolumn")
parser.add_argument("-u", action="store", dest="unaligned")
parser.add_argument("-con", action="store", dest="conversion")
parser.add_argument("filestoalign", nargs=argparse.REMAINDER)
arguments = parser.parse_args()

output = []
refcolumns = arguments.refcolumns.split(",")
rac = arguments.refalignmentcolumn

if arguments.conversion:
	with open(arguments.conversion) as file:
		r = [word.strip() for word in file]
		r2 = ( line.split('\t') for line in r )
		dict = { row[0]:row[1:] for row in r2 }

for refrow in range(len(refcolumns)):
	if refcolumns[refrow] == rac:
		refalignmentcolumn = refrow

with open(arguments.referencefile, 'rU') as csvfile:
	reflines = csv.reader(csvfile)
	for refline in reflines:
		temp = []
		for refcolumn in refcolumns:
			temp.append(refline[int(refcolumn)])
		output.append(temp)

length = len(output)
position = len(refcolumns) + 1

for argument in arguments.filestoalign:
	file = argument.split(":")[0]
	columns = argument.split(":")[1].split(",")
	ac = argument.split(":")[2]
	for row in range(len(columns)):
		if columns[row] == ac:
			alignmentcolumn = row
	
	list = []
	with open(file, 'rU') as csvfile:
		lines = csv.reader(csvfile)
		for line in lines:
			temp = []
			for column in columns:
				temp.append(line[int(column)])
			list.append(temp)	
			
	if arguments.conversion:
		header = 0
		
		for i in list:
			if header == 0:
				for blank in range (position - len(output[0])):
					output[0].append("")
				for head in i:
					output[0].append(head)
				header = 1
				
			else:
				matched = 'false'
				for r in output[1:length]:
					checking = dict.get(r[refalignmentcolumn], None)
					if checking is not None:
						if i[alignmentcolumn] in checking:
							for vi in range (position - len(r)):
								r.append("")
							for vii in i:
								r.append(vii)
							matched = 'true'
				
				if matched == 'false':
					if arguments.unaligned:
						if len(output) == length:
							output.append([])
						if len(output[len(output) - 1]) == position + len(columns):
							output.append([])
							for ii in range (position):
								output[len(output) - 1].append("")
							for iii in i:
								output[len(output) - 1].append(iii)
						
						else:
							for unalrow in range(length,len(output)):
								if (len(output[unalrow]) == position + len(columns)):
									continue
								else:
									for iv in range (position - len(output[unalrow])):
										output[unalrow].append("")
									for v in i:
										output[unalrow].append(v)
									break
	
	else:
		
		header = 0
		
		for i in list:
			if header == 0:
				for blank in range (position - len(output[0])):
					output[0].append("")
				for head in i:
					output[0].append(head)
				header = 1
				
			else:
				matched = 'false'
				for r in output[1:length]:
					if i[alignmentcolumn] == r[refalignmentcolumn]:
						for vi in range (position - len(r)):
							r.append("")
						for vii in i:
							r.append(vii)
						matched = 'true'
				
				if matched == 'false':
					if arguments.unaligned:
						if len(output) == length:
							output.append([])
						if len(output[len(output) - 1]) == position + len(columns):
							output.append([])
							for ii in range (position):
								output[len(output) - 1].append("")
							for iii in i:
								output[len(output) - 1].append(iii)
							
						else:
							for unalrow in range(length,len(output)):
								if (len(output[unalrow]) == position + len(columns)):
									continue
								else:
									for iv in range (position - len(output[unalrow])):
										output[unalrow].append("")
									for v in i:
										output[unalrow].append(v)
									break					
	
	position += len(columns) + 1
	
with open(arguments.output, "wb") as csvfile:
	writer = csv.writer(csvfile)
	writer.writerows(output)
	
#test: python alignment.py -r 19F_SDMM_PAUL.csv -o testing.csv -u 1 -c 2,7,3,8 -ac 2 19F_SDMM_PIPE1.csv:0,1,6,3:0 19F_SDMM_PIPE2.csv:0,1,6,3:0
#test: python alignment.py -r 19F_SDMM_PAUL.csv -o testing2.csv -u 1 -c 2,7,3,8 -ac 2 -con TIGR4vsD39test.dat D39_SDMM_T2.csv:0,1,6,3:0