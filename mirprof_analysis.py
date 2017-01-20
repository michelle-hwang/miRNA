from __future__ import print_function
from Bio import SeqIO
import argparse
import sys
from operator import itemgetter

# -----------------------------------------------------------------------------

def get_counts(ranking, args):
	"""Gets count data for miRNAs whose origin species is in the priority list.
	Args:
		ranking: Dict in which keys = rank #, values = species
	Returns:
		A dict in which keys = miRNA and values are lists. Lists 
		include species and count data. 

		EXAMPLE:
		{'miR10': [[['species1', 10, 5.8, 2.3], ['species2', 54, 54, 15.5]}
	"""

	try:
		infile = open(args.infile, 'r')
	except:
		print('ERROR: Cannot open ', args.infile, '!', sep='')

	db = dict()
	species = 0

	for line in infile:
		if "Organism:" in line:
			species = line.split()[1][:-1]
		elif ("miR" in line) and ("weighted" not in line) and (species in ranking.values()):
			cols = line.rstrip()
			cols = cols.split(",")
			mirna = cols[0].split('"')[1]
			if mirna in db:
				db[mirna].append([species, cols[1], cols[2], cols[3]])
			else:
				db[mirna] = [species, cols[1], cols[2], cols[3]]

	infile.close()
	print()
	return db


def get_seq(name, args):
	"""Given miRNA ID, gets FASTA sequence.

	Uses SeqIO in BioPython to parse the FASTA infile. Looks for a 
	match in name to FASTA header.

	Args:
		name: miRNA ID
	Returns:
		string: FASTA sequence
	"""

	with open(args.fasta) as fh:
		for record in SeqIO.parse(fh, 'fasta'):
			if record.id in name or name in record.id:
				return(record.seq)


def initialize_priority_lists(plist):
	"""Obtains species priority information.
	Args: 
		prioritylist: file handle of input prioirty list 
	Returns:
		dict: key = rank #, value = species
	"""

	try:
		priority_list = open(plist, 'r')
	except:
		print("ERROR: Cannot open ", plist, "!", sep="")

	ranking = dict()
	num_rank = 1 # Key in dict = rank
	for species in priority_list:
		s = species.rstrip()
		ranking[num_rank] = s.lower()
		num_rank += 1
	priority_list.close()
	return ranking


def pick_best_mirna(candidates, ranking):
	"""Chooses a species as top candidate for an miRNA based on priority list.

	Loops over ranked species by descending order (1st, 2nd...) and searches
	the candidate list for a matching species. If there is a match, returns 
	tuple in candidate list.

	Args:
		candidates: list of tuples with species and count information
		ranking: dict of species with priority
	Returns:
		tuple: top species with count information
	"""
	for x in range(1, len(ranking), 1):
		if(len(candidates) == 4):
			return candidates
		for a in candidates:
			if ranking[x] == a[0]:
				return a


def print_mirna(key, info, args, outfile):
	"""Prints final miRNAs in a tab-delimited .txt file."""
	print(key, end='\t')
	print('\t'.join([x for x in info]), get_seq(key, args), sep="\t")

	print(key, end='\t', file=outfile)
	print('\t'.join([x for x in info]), get_seq(key, args), sep="\t", file=outfile)


# -----------------------------------------------------------------------------

def main():

	# Command line parameters and help
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
		description='''
	This script will take mirProf output in which species are not 
	grouped/collapsed and a ranked list of species in which priority is 
	given. A representative sequence will be determined for each miRNA 
	based on the priority list. 

	Default usage:

		prioritize_mirprof.py sequences.fa myresults.csv priority_list.txt outfile_name

	///AUTHOR: Michelle Hwang
	///DATE: 12/11/2016''')
	parser.add_argument('fasta', help = 'Name of fasta file output from mirProf.')
	parser.add_argument('infile', help = 'Name of .csv file output from mirProf.')
	parser.add_argument('plist', help = '''Name of text file containing 
		priority list for species. One species per line. First species has 
		highest priority. Please use three letter code for naming species as 
		determined by mirbase.''')
	parser.add_argument('outfile', help = '''Specify an output file name.''')
	args = parser.parse_args()

	# Open and print to outfile
	try:
		args.outfile = args.outfile.strip()
		outfile = open(args.outfile, 'w')
	except:
		print('ERROR: Cannont open ', args.outfile, '!', sep="")

	# Get species priority data from files
	print('\nBEGIN: Importing species from priority list.')
	rank = initialize_priority_lists(args.plist)
	for k, v in rank.iteritems():
		print(k, v, sep="\t")

	# Get miRNA count information from .csv file
	print('\nBEGIN: Importing miRNA count information from .csv file.')

	db = get_counts(rank, args)
	if not db:
		print('ERROR: miRNA count db empty. Process halted.')
		raise SystemExit
	else:
		print('Successful!')

	# Get best miRNA hit based on priority rank
	print('\nBEGIN: Determining best miRNA species candidates.')
	for key, value in db.iteritems():
		best_mirna_info = pick_best_mirna(value, rank)
		print_mirna(key, best_mirna_info, args, outfile)

	outfile.close()


main()
