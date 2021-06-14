#!/usr/bin/env python
import sys, gzip, os
import getopt
from glob import glob
from time import sleep
import re

# This script removes the first 8bp of a read and adds the sequence to the readID. The quality information is discarded

# This is the expected structure of the FastQ files:

# barcode (UMI) (8bp)    //    PolyT     //      Insert

# After moving the UMI sequences, the script looks for up to 3 T at the start of the sequence, and removes those.
# Sequences with more than 3 Ts at the 5' end are clipped a maximum of 3 TTT

# Script last modified 21 July 2020, Felix Krueger

polyT = {}         # storing the number of Poly Ts at the start of the read (after the UMI)
fhs = {}           # storing the filehandles for all output files

def submain():
	
	print (f"Python version: {sys.version}.")
	allfiles = glob("*.fastq.gz")
	# print (allfiles)
	allfiles.sort() # required as glob doesn't necessarily store files in alphabetical order
	# print (allfiles)

	for filename in allfiles:
		print (f"Reading in FastQ file:\t >> {filename} <<\n")
		main(filename)
		polyT = {} # resetting

def main(filename):

	expected_count   = 0
	unexpected_count = 0   
	count = 0              # total sequence count
	
	print (f"Reading file: >{filename}<")
	
	outfh = make_out_filehandle("UMIed",filename)
	pattern = '^T+'
	p = re.compile(pattern)

	with gzip.open(filename) as cf:
	
		while True:
			readID  = cf.readline().decode().strip()
			seq     = cf.readline().decode().strip()
			line3   = cf.readline().decode().strip()
			qual    = cf.readline().decode().strip()
			
			if not qual:
				break
			
			count += 1

			if count%500000 == 0:
				print (f"Processed {count} reads so far")

			# print (f"{readID}\n{seq}\n{line3}\n{qual}\n")
		
			## STEP 1: Remove UMI and write it into the read ID

			# These are the expected in-line codes:
			# barcode (8bp)    //    PolyT     //      Insert
		
			barcode   = seq[0:8]
			rest      = seq[8::]
			qual_rest = qual[8::]
			#print (f"sequence: {seq}\nbarcode:  {barcode}\nrest:             {rest}\n")
			readID += f':{barcode}'
			#print (f"{readID}\n{rest}\n{line3}\n{qual_rest}\n")
			#sleep(1)


			## STEP 2: Now we need to find a number of PolyTs at the start

			m = p.match(rest)

			
			if m is None:
				# Using what we aleady have if the read does not have T(s) at the start
				new_rest = rest
				new_rest_qual = qual_rest
				
			else:
				polyTlength = len(m[0])
				# print (m[0])
				if not m[0] in polyT: # This is to keep track of the number of T(s) trimmed
					polyT[m[0]] = 0

				# removing only up to 3 bp of T to avoid genomic polyA trimming
				if polyTlength > 3:
					polyTlength = 3
					if not 'TTT' in polyT:
						polyT['TTT'] = 0
						
					polyT['TTT'] += 1  
				else:
					polyT[m[0]] += 1


				# removing polyT and its quality scores
				new_rest      = rest[polyTlength::]
				new_rest_qual = qual_rest[polyTlength::]
				
			readID = readID.replace(" ","_") # this is required for e.g. Bowtie2 to retain the last part of the read ID (= the UMI sequence)

			
			# print ("\n".join([readID, new_rest, line3, new_rest_qual]))
			outfh.write (("\n".join([readID, new_rest, line3, new_rest_qual]) + "\n").encode())
	
	outfh.close()
	
	print (f"Total number of reads processed: {count}")
	for tees in sorted (polyT.keys()):
		print (f"{tees}\t{polyT[tees]}")

def make_out_filehandle(sample_name,filename):
	
	print (f"Got following sample name: {sample_name}")
	
	# extracting useful parts from filename
	# Example name: lane7265_ACTTGA_fob1_YPD_LIGseq_L001_R1.fastq.gz

	pattern = '(lane.*_L00\d)_(R\d.fastq.gz)'
	p = re.compile(pattern)
	print (filename)
	m = p.findall(filename)
	sample = m[0][0]
	ending = m[0][1]
	new_filename = f"{sample}_{sample_name}_{ending}"
	# print (new_filename)
	
	outfh  = gzip.open (new_filename,mode='w',compresslevel=3)
	
	return outfh
	
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)	


if __name__ == "__main__":
	submain()
else:
	print ("Just getting imported")
