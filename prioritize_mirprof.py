
import argparse
import sys
from operator import itemgetter


## Command line parameters and help
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
	description='''
This script will take mirProf output in which species are not 
grouped/collapsed and a ranked list of species in which priority is 
given. A representative sequence will be determined for each miRNA 
based on the priority list. The output fasta file will contain only
the representative sequence for each miRNA with the following header:

	>miRNA-name_species-raw_count-normalized_weighted_count

Default usage:

	python3 prioritize_mirprof.py sequences.fa myresults.csv plant_list.txt animal_list.txt outfile_name

///AUTHOR: Michelle Hwang
///DATE: 10/6/2015''')
parser.add_argument('fasta', help = 'Name of fasta file output from mirProf.')
parser.add_argument('infile', help = 'Name of .csv file output from mirProf.')
parser.add_argument('priority_list_p', help = '''Name of text file containing 
	priority list for PLANT species. One species per line. First species has 
	highest priority. Please use three letter code for naming species as 
	determined by mirbase.''')
parser.add_argument('priority_list_a', help='''Name of text file containing
	priority list for ANIMAL speces. One species per line. First species has
	highest priority. Please use three letter code for naming species as
	determined by mirbase.''')
parser.add_argument('outfile', help = '''Specify an output file name.''')
parser.add_argument('-nm', '--no_match', type=int, default=1, help='''If miRNA
	does not have a match to any of the species on the priority list: 
	(1) Discard = DEFAULT
	(2) Choose sequence with most raw counts''')
parser.add_argument('-a', '--alternative_outfile', type=bool, default=False,
	help='''Instead of a fasta file output, the output will consist of a tab
	delimited table. Turn on by setting to "True"''')


args = parser.parse_args()
fasta = open(args.fasta, 'r')
infile = open(args.infile, 'r')
plant_list = open(args.priority_list_p, 'r')
animal_list = open(args.priority_list_a, 'r')
outfile = open(args.outfile, 'w')
no_match = args.no_match
alt_out = args.alternative_outfile


## Get PLANT priority list
ranking_p = dict()
num_rank = 1 # Key in dict = rank
for species in plant_list:
	ranking_p[num_rank] = species.rstrip()
	num_rank += 1
plant_list.close()
print("Plant priority list:")
for k,v in ranking_p.items():
	print(k,v)


## Get ANIMAL priority list
ranking_a = dict()
num_rank = 1
for species in animal_list:
	ranking_a[num_rank] = species.rstrip()
	num_rank += 1
animal_list.close()
print("", "Animal priority list:", sep="\n")
for k,v in ranking_a.items():
	print(k,v)


db = {} # Dict of {miRNA, [(species, seq, count)]
mirna_db = [] # List of tuple: (species, seq, count)


## Populate the db; go through entire file
for line in fasta:
	if( line.startswith('>')):
		cols = line[1:].split('_')
		species = cols[0][:3]
		mirna = cols[0][4:]
		count = cols[2][0]
	else:
		if(mirna in db):
			db[mirna].append((species, count, line.rstrip()))
		else:
			db[mirna] = [(species, count, line.rstrip())]
fasta.close()


def get_counts(name, species):
	hit_species = False
	for line in infile:
		if species in line:
			hit_species = True
		if name in line and hit_species is True:
			values = line.split(',')
			infile.seek(0,0)
			return(values[1], values[3].rstrip()) # Return raw count and weighted, normalized count
	print("ERROR: Hit end of .csv file without finding any miRNA hits.")


## Print a header if tab-delimited output
if(alt_out is True):
	print("mirna", "species", "raw_count", "norm_count", "seq", sep='\t', end='\n', file=outfile)


## Pick the best miRNA sequence based on priority
for mirna, mirna_list in sorted(db.items()):

	is_ranked = False
	## If animal mirna:
	if '-' in mirna: 
		for x in range(1,len(ranking_a)): # Look through priority list 
			matches = [items for items in mirna_list if ranking_a[x] == items[0]]
			if matches: 
				top_match = max(matches,key=itemgetter(1)) # Use sequence with most raw counts
				raw_count, norm_count = get_counts(mirna, top_match[0])
				if(alt_out is True):
					print(mirna, top_match[0], raw_count, norm_count, top_match[2], sep="\t", end='\n', file=outfile)
				else:
					print('>'+mirna+'_'+top_match[0]+'_'+raw_count+'_'+norm_count, top_match[2], sep="\n", end="\n", file=outfile)
				is_ranked = True
				break

	## If plant mirna:
	else: 
		for x in range(1,len(ranking_p)): 
			matches = [items for items in mirna_list if ranking_p[x] == items[0]]
			if matches:
				top_match = max(matches,key=itemgetter(1))
				raw_count, norm_count = get_counts(mirna, top_match[0])
				if(alt_out is True):
					print(mirna, top_match[0], raw_count, norm_count, top_match[2], sep="\t", end='\n', file=outfile)
				else:
					print('>'+mirna+'_'+top_match[0]+'_'+raw_count+'_'+norm_count, top_match[2], sep="\n", end="\n", file=outfile)
				is_ranked = True
				break

	## If mirna has no matches to any priority list and you do not want to discard:
	if(is_ranked == False and no_match is 2): 
		top_match = max(mirna_list,key=itemgetter(1))
		raw_count, norm_count = get_counts(mirna, top_match[0])
		if(alt_out is True):
			print(mirna, top_match[0], raw_count, norm_count, top_match[2], sep="\t", end='\n', file=outfile)
		else:
			print('>'+mirna+'_'+top_match[0]+'_'+raw_count+'_'+norm_count, top_match[2], sep="\n", end="\n", file=outfile)


infile.close()
outfile.close()


