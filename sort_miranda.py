#!/usr/local/python/3.4.0/bin/python3

import argparse
import os
import operator
import re
import subprocess as sub
import time

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description='''
This script will take in grep '>>' miranda output and 
filter out and sort sequences based on:
	1. Energy threshold
	2. Score threshold
	3. Near the 3'UTR 

- Assumes there is only one location per transcript target

Default usage:

	sort_miranda.py output.miranda transcripts.fa outfile.txt

	TEMP: python3 ~/scripts/sort_miranda.py ctl.mirprof.fa.miranda.short YSA_transdecoder.gff3 all.miranda.sorted
	TEST: python3 ~/scripts/sort_miranda.py test.miranda YSA_transdecoder.gff3 test.miranda.out -s 150

///AUTHOR: Michelle Hwang
///DATE: 7/6/2016''')
parser.add_argument('miranda', help = 'Name of output from miranda.')
parser.add_argument('transdecoder', help = 'Name of transdecoder output of transcriptome')
parser.add_argument('outfile', help = 'Name for output file')
parser.add_argument('-e', '--energy', type=int, default=20, help='''
	Energy threshold (default=20)
	- e.g. remove those folds with larger energy values than |20|
	- Most likely secondary structure has min free energy''')
parser.add_argument('-s', '--score', type=int, default=100, help='''
	Score threshold (default=100)
	- e.g. remove those matches with scores less than 100''')

args = parser.parse_args()
wkdir = os.getcwd()

# ------------------------------------------------------------------------------------------------

class MyError(Exception):
	pass

class Mirna:
	def __init__(self, name):
		self.name 			= name
		self.coordinates 	= list()
		self.energies 		= list()
		self.scores 		= list()
		self.targets 		= list()
		self.ranks 			= list()

	def add_target(self, target, energy, score, coordinate):
		self.targets.append(target)
		self.energies.append(energy)
		self.scores.append(score)
		self.coordinates.append(coordinate)

	def rank_targets(self):
		self.ranks = [operator.itemgetter(0)(x) for x in sorted(enumerate(self.scores,1), key=operator.itemgetter(1))]

	def print_mirna(self):
		print(self.name)


# ------------------------------------------------------------------------------------------------


def print_time(start):
	"""Prints time in hours, minutes, seconds since program elasped."""
	stop = time.clock()
	m, s = divmod(stop-start, 60)
	h, m = divmod(m, 60)
	print( '%d:%02d:%02d' % (h, m, s))


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
		print_time(start_time)
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

	print(">>> 3'UTR information collected.")
	print_time(start_time)
	return threeprimes


def check_score(score):
	"""Returns false if target does not pass score threshold.
	True otherwise."""
	if args.score > score:
		return False
	return True


def check_energy(energy):
	"""Returns false if target does not pass energy threshold.
	True otherwise."""
	if args.energy < abs(energy):
		return False
	return True


def check_coordinate(coordinate, threeprime, transcript_len):
	"""Returns false if target does is not in 3'UTR.
	Take in coordinate (length, pos) and threeprime info
	in form of (start, end, strand, length).
	True otherwise."""

	target_start = threeprime[0]
	target_end	 = threeprime[1]
	strand 		 = threeprime[2]

	for hit in coordinate:
		if strand is '+':
			mirna_start = hit[0]
			mirna_end 	= hit[1]
		elif strand is '-':
			mirna_start = transcript_len - mirna_end
			mirna_end 	= transcript_len - mirna_start
		else:
			raise MyException('Invalid strand. Must be + or -.')		

		if ((mirna_start >= target_start) and (mirna_end <= target_end)):
			return True
	
	print('\tWARNING: No hits were found in the 3 prime UTR region.')
	return False


def print_outfile(mirnas):
	"""Formats and prints result outfile."""
	outfile = open(args.outfile, 'w')
	for m in mirnas: # for each mirna
		for n in range(1,len(mirnas[m].targets),1):
			co = mirnas[m].coordinates[n]
			print(m, mirnas[m].ranks[n],
					 mirnas[m].targets[n], 
					 ", ".join([str(i[0]) for i in co]),
					 mirnas[m].energies[n],
					 mirnas[m].scores[n],
					 sep="\t", end="\n", file=outfile)
	outfile.close()

# ------------------------------------------------------------------------------------------------

def main():

	all_mirnas		= dict() # Mirna name, class object, STATIC
	current_mirna 	= None 
	data			= None # Temp container of Mirna object  

	threeprimes = get_transdecoder_info()
	miranda 	= open(args.miranda, 'r')
	lines 		= miranda.readlines()

	miranda.close()
	print(">>> Miranda output collected.")
	print_time(start_time)

	print(">>> Reading Miranda lines.")
	for line in lines:
		line 	= line.rstrip()
		fields 	= line.split('\t')

		mirna 			= fields[0]
		mirna			= mirna[2:] # Remove '>>'
		target 			= fields[1]
		score 			= float(fields[2])
		energy 			= float(fields[3])
		length 			= int(fields[7]) # Length of mirna
		transcript_len  = int(fields[8]) # Length of target transcript
		pos 		    = fields[9].lstrip()
		pos				= pos.split() # Can have more than one
		pos 			= list(map(int, pos)) # Change to int

		coordinate = list()
		for p in pos:
			coordinate.append((p, p+length))

		# If target does not pass energy or score threshold
		if not check_score(score):
			continue
		if not check_energy(energy):
			continue

		# If target location is not in the 3' UTR region
		for c in coordinate:
			if ((c in threeprimes) and
				(check_coordinate(c, threeprimes[target], transcript_len) is False)):
				coordinate.remove(c)
		if coordinate is None:
			continue

		# Build Mirna object
		if current_mirna is None:
			current_mirna = mirna
			data = Mirna(current_mirna)
		elif current_mirna != mirna:
			data.rank_targets()
			all_mirnas[data.name] = data # Add Mirna object
			current_mirna = mirna
			data = Mirna(current_mirna)
		else:
			pass
		data.add_target(target, energy, score, coordinate)

	data.rank_targets()
	all_mirnas[data.name] = data

	print(">>> Target filtering and ranking completed.")
	print_outfile(all_mirnas)


if __name__ == "__main__":
	start_time = time.clock()
	main()
	print(">>> Script complete.")
	print_time(start_time)







