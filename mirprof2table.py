from __future__ import print_function
from __future__ import division
import argparse
import sys
# Requires Python 3 or newer


parser = argparse.ArgumentParser(description='''Takes in FASTA file output 
	from mirProf and converts it to tab-delimited table format.''')
parser.add_argument('infile', help='Name of fasta file output from mirProf.')
parser.add_argument('outfile', help='''Specify an output file name.''')
parser.add_argument('-g', '--grouped', type=bool, default=False, help='''Were all oragnisms
	grouped in the outfile? True or False.''')
parser.add_argument('-l', '--length', type=int, default=0, help='''Remove 
	any sequences that have a raw count of this length or lower. Rule does not apply
	if collapse is true. Default=0''')
parser.add_argument('-c', '--collapse', type=bool, default=False, help='''True or False. 
	Collapse miRNA with the same match, despite different sequence. Will not collapse
	organisms. Default=False''')
parser.add_argument('-al', '--add_length', type=bool, default=False, help='''True or
	False. Adds an additional column that specifies sequence length. Does not apply
	if collapse is true. Default=False''')

args=parser.parse_args()
infile = open(args.infile, 'r')
outfile = open(args.outfile, 'w')
grouped_or_not = args.grouped
l = int(args.length)
collapse = args.collapse
add_length = args.add_length


current_mirna = None
current_count = 0

print("grouped_or_not: ", grouped_or_not)
print("collapse: ", collapse)
print("add_length: ", add_length)

#Print header
if(grouped_or_not is True):
	if(collapse is True):
		print("mirna", "count", sep="\t", end="\n", file=outfile)
	else:
		if(add_length is True):
			print("mirna", "consecutive_match", "count", "length", "sequence", sep="\t", end="\n", file=outfile)
		else:
			print("mirna", "consecutive_match", "count", "sequence", sep="\t", end="\n", file=outfile)
else:
	if(collapse is True):
			print("species", "mirna", "count", sep="\t", end="\n", file=outfile)
	else:
		if(add_length is True):
			print("species", "mirna", "consecutive_match", "count", "length", "sequence", sep="\t", end="\n", file=outfile)
		else:
			print("species", "mirna", "consecutive_match", "count", "sequence", sep="\t", end="\n", file=outfile)




for line in infile:

	if(line.startswith('>')):
		if( grouped_or_not is True ): 

			#GROUPED = TRUE
			#EX: >all_combined-let-7_1_1x

			cols = line[14:].split('_')

			if(collapse is True):				
				if(current_mirna is None): 
					#If first mirna
					current_mirna = cols[0]
					current_count = int(cols[2][0])

				elif(str(cols[0]) != str(current_mirna)): 
					#If new mirna
					print(current_mirna, current_count, sep="\t", end="\n", file=outfile)
					current_mirna = cols[0]
					current_count = int(cols[2][0])

				else: 
					#If same mirna
					current_count = current_count + int(cols[2][0])
			else:
				if(l < int(cols[2][0])):
					print(cols[0], cols[1], cols[2][0], sep="\t", end="\t", file=outfile)
				else:
					continue

		else: 

			#GROUPED = FALSE
			#EX: >aae-let-7_1_1x

			cols = line[1:].split('_')
			species = cols[0][:3]
			name = cols[0][4:]

			if(collapse is True):
				if(current_mirna is None): 
					#If first mirna
					current_mirna = name
					current_count = int(cols[2][0])

				elif(str(name) != str(current_mirna)): 
					#If new mirna
					print(species, current_mirna, current_count, sep="\t", end="\n", file=outfile)
					current_mirna = name
					current_count = int(cols[2][0])

				else: 
					#If same mirna
					current_count = current_count + int(cols[2][0])				

			else:
				if( l < int(cols[2][0])):	
					print(species, name, sep="\t", end="\t", file=outfile)
					print(cols[1], cols[2][0], sep="\t", end="\t", file=outfile)
				else:
					continue

	else:
		# Will not print sequence if collapse is true
		if(collapse is False and l < int(cols[2][0])):
			if( add_length==True):
				print(len(line.strip()), line.strip(), sep="\t", end="\n", file=outfile)
			else:
				print(line.strip(), end="\n", file=outfile)

#Gets skipped last sequence
if(collapse is True and grouped_or_not == 1):				
	print(current_mirna, current_count, sep="\t", end="\n", file=outfile)
if(collapse is True and grouped_or_not == 2):
	print(species, current_mirna, current_count, sep="\t", end="\n", file=outfile)



infile.close()
outfile.close()



