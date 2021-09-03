#!/usr/bin/env python
import sys, gzip, os
import getopt
from glob import glob
from time import sleep
import re

# This script adds the first 23bp of a read sequence to the readID. 
# It then uses these 23bp, as an in-memory sequence-based deduplication step. 
# This is an experimental pre-preprocessing step to evaluate the influence
# of multi-mapping loci in TrAEL-seq data 

# Script last modified 20 Jan 2021, Felix Krueger

## fhs = {}           # storing the filehandles for all output files

def submain():
	
	print (f"Python version: {sys.version}.")
	# allfiles = glob("*.fastq.gz") 
	allfiles = sys.argv[1:] # enables loop processing
	allfiles.sort() # required as glob doesn't necessarily store files in alphabetical order

	for filename in allfiles:
		print (f"Reading in FastQ file:\t >> {filename} <<\n")
		main(filename)
		polyT = {} # resetting

def main(filename):

	expected_count   = 0
	unexpected_count = 0   
	count = 0              # total sequence count
	unique_count = 0
	dups_count = 0
	observed = {}  # dict storing all sequences (23bp) we have seen already

	print (f"Reading file: >{filename}<")
	
	outfh = make_out_filehandle("seqDedup",filename)
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
		
			barcode   = seq[0:23]
			# print (f"sequence: {seq}\nbarcode:  {barcode}\n")
			readID += f':{barcode}'
			# print (f"{readID}\n{seq}\n{line3}\n{qual}\n")


			if barcode in observed.keys():
				# print ("already present, skipping...")
				dups_count += 1
				continue
			else:
				observed[barcode] = 1
				unique_count += 1
				readID = readID.replace(" ","_") # this is required for e.g. Bowtie2 to retain the last part of the read ID (= the UMI sequence)
				# print ("\n".join([readID, new_rest, line3, new_rest_qual]))
				outfh.write (("\n".join([readID, seq, line3, qual]) + "\n").encode())

			
	
	outfh.close()
	
	print (f"Sequences processed: {count}\nUnique sequences: {unique_count}\nDuplicate sequences (removed): {dups_count}\n")

	
def make_out_filehandle(sample_name,filename):
	
	print (f"Got following sample name: {sample_name}")
	
	# extracting useful parts from filename
	# Example name: lane7265_ACTTGA_fob1_YPD_LIGseq_L001_R1.fastq.gz
	# lane7360_CAGATC_Colo_2_day_SEL_L001_R1.fastq.gz
	pattern = '(lane.*)_(L00\d)_(R\d.fastq.gz)'
	p = re.compile(pattern)
	print (filename)
	m = p.findall(filename)
	sample = m[0][0]
	lanenumber = m[0][1]
	ending = m[0][2]
	new_filename = f"{sample}_{sample_name}_{lanenumber}_{ending}"
	# print (new_filename)
	
	outfh  = gzip.open (new_filename,mode='wb',compresslevel=3)
	
	return outfh
	
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)	


if __name__ == "__main__":
	submain()
else:
	print ("Just getting imported")
