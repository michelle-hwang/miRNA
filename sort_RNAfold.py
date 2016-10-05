#!/usr/local/python/3.4.0/bin/python3

import subprocess as sub
import argparse
import re
import os
import operator

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description='''
# Currently runs for 30 min for 6000 hits

This script will take in a target-predict coordinates for predicted target sites of
miRNA and RNAfold output of predicted secondary structures of the target mRNAs. It
will determine whether the binding site of the miRNA on the target mRNA contains
a loop or whether it is a stem.

Default usage:

	python3 sort_RNAfold.py rnafold.fa coordinates.csv transdecoder.gff
	python3 ~/scripts/sort_RNAfold.py Mxg.targets.rnafold Mxg.targets.csv MXG_transdecoder.gff3 > test

///AUTHOR: Michelle Hwang
///DATE: 6/29/2016''')
parser.add_argument('fasta', help 	= 'Name of hairpin file output from RNAfold.')
parser.add_argument('infile', help 	= 'Name of .csv file output from target-predict.')
parser.add_argument('transdecoder', help = 'Name of .gff Transdecoder output.')

args 	= parser.parse_args()
infile 	= open(args.infile, 'r')
fasta 	= open(args.fasta, 'r')

hairpins = fasta.readlines()
fasta.close()	
print(">>> RNAfold output information collected.")
wkdir = os.getcwd()

# ------------------------------------------------------------------------------------------------

class Mirna:
	def __init__(self, name):
		self.name = name
		self.coordinates = list()
		self.energies = list()
		self.structures = list()
		self.targets = list()
		self.ranks = list()

	def add_target(self, target, coordinate):
		self.targets.append(target)
		self.coordinates.append(coordinate)

	def determine_folds(self):
		for i in range(0,len(self.targets),1):
			(structure, energy) = determine_fold(self.targets[i], self.coordinates[i])
			self.structures.append(structure)
			self.energies.append(energy)

	def rank_targets(self):
		self.ranks = [operator.itemgetter(0)(x) for x in sorted(enumerate(self.structures,1), key=operator.itemgetter(1))]

	def print_mirna(self):
		print(self.name)
		
# ------------------------------------------------------------------------------------------------

def get_transdecoder_info():
	"""Takes transdecoder .gff3 output and takes 3'UTR information.
	Returns a dictionary of transcripts with values in the form of
	a tuple (pos1, pos2, strand + or -). Will only take best 3'UTR 
	region."""

	print(">>> Attempting to gather Transdecoder .gff3 info using grep.")

	try:
		lines = sub.getoutput("grep 'three_prime_UTR' "+args.transdecoder)
		lines = lines.split('\n')
		print(">>> Grep from Transdecoder gff file successful.")
	except IOError:
		print('\tERROR: Transdecoder .gff3 file could not be found.')

	threeprimes = dict()
	duplicates = list()
	current_transcript = None

	for line in lines:
		line 	= line.rstrip()
		fields 	= line.split('\t')

		transcript 	= fields[0]
		pos1		= int(fields[3])
		pos2 		= int(fields[4])
		strand		= fields[6]

		# Remove unlikely 3'UTRs
		if (pos2-pos1) < 25:
			continue 
		elif current_transcript is transcript:
			duplicates.append((pos1, pos2, strand, pos2-pos1))
		else:
			if current_transcript is not None:
				# Picks hit with shortest 3'UTR
				top_hit = sorted_dups = sorted(duplicates, key=operator.itemgetter(0))[0]
				threeprimes[transcript] = top_hit
				duplicates = list()
			current_transcript = transcript 
			duplicates.append((pos1, pos2, strand, pos2-pos1))

	return threeprimes

def check_coordinate(coordinate, threeprime, transcript_len):
	"""Returns false if target does is not in 3'UTR.
	Take in coordinate (length, pos) and threeprime info
	in form of (start, end, strand, length).
	True otherwise."""

	target_start = int(threeprime[0])
	target_end	 = int(threeprime[1])
	strand 		 = threeprime[2]
	coordinates	 = coordinate.split('-')
	mirna_start	 = int(coordinates[0])
	mirna_end	 = int(coordinates[1])

	if strand is '+':
		pass
	elif strand is '-':
		mirna_start = transcript_len - mirna_end
		mirna_end 	= transcript_len - mirna_start
	else:
		raise MyException('Invalid strand. Must be + or -.')		

	if ((mirna_start >= target_start) and (mirna_end <= target_end)):
		return True
	
	return False


def determine_fold(target, coordinates):
	dotbracket 		= get_dotbracket(target)
	(start, stop) 	= coordinates.split('-')

	d = dotbracket[0] #dotbracket
	e = dotbracket[1] #energy
	e = e.strip()
	e = re.sub('[()]', '', e)
	s = d[int(start)-1:int(stop)-1] #subseq

	if re.findall(r'\(+|\)+',s):
		max_consec_loops = max(len(n) for n in re.findall(r'\(+|\)+',s))
	else: 
		max_consec_loops = 0

	if re.findall(r'(.)+',s):
		max_consec_stems = max(len(n) for n in re.findall(r'(.)+',s))
	else:
		max_consec_stems = 0

	if(s.count('.') is 0):
		ratio = s.count('(')+s.count(')')
	else:
		ratio = (s.count('(')+s.count(')'))/s.count('.')

	structure = (max_consec_loops, max_consec_stems, ratio)
	return(structure, e)

def get_dotbracket(target):
	out = sub.check_output("awk '/"+target+"$/ {print FNR}' "+wkdir+"/"+args.fasta+" | head -n 1", universal_newlines=True, shell=True)
	return(hairpins[int(out)+1].split(' '))

def get_length(target):
	out = sub.check_output("awk '/"+target+"$/ {print FNR}' "+wkdir+"/"+args.fasta+" | head -n 1", universal_newlines=True, shell=True)
	return(len(hairpins[int(out)]))

def print_out(all_mirnas):
	for m in all_mirnas: # for each mirna
		for n in range(1,len(all_mirnas[m].energies),1):
			print(m, all_mirnas[m].ranks[n],
					 all_mirnas[m].targets[n], 
					 all_mirnas[m].coordinates[n], 
					 all_mirnas[m].energies[n],
					 all_mirnas[m].structures[n][0],
					 all_mirnas[m].structures[n][1],
					 all_mirnas[m].structures[n][2],
					 sep="\t", end="\n")

# ------------------------------------------------------------------------------------------------

def main():

	all_mirnas 		= dict() # all mirna names and their class objects
	current_mirna 	= None 
	threeprimes		= get_transdecoder_info()
	print(">>> Transdecoder information collected.")

	for line in infile:

		if "No target found" in line or "sRNA ID" in line:
			continue 

		spline 		= line.split(',')
		mirna		= spline[0]
		target 		= spline[1]
		coordinate 	= spline[2]

		# If target location is not in the 3'UTR region:
		if ((target in threeprimes) and (check_coordinate(coordinate, threeprimes[target], get_length(target)) is False) or coordinate is None):
			continue	

		if(spline[0] == current_mirna):
			data.add_target(target, coordinate)

		else:
			if(current_mirna is not None):
				data.determine_folds()
				data.rank_targets()
				all_mirnas[data.name] = data
			current_mirna = mirna
			data = Mirna(current_mirna)
			data.add_target(target, coordinate)

	current_mirna = mirna
	data = Mirna(current_mirna)
	data.add_target(target, coordinate)
	infile.close()

	print_out(all_mirnas)


if __name__ == "__main__":
    main()

