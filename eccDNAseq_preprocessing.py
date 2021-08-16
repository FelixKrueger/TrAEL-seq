#!/usr/bin/env python
import sys, gzip, os
import getopt
from glob import glob
from time import sleep
import re

# This script removes the first 8bp of a read and adds the sequence to the readID. The quality information is discarded

# This is the expected structure of the FastQ files:

# barcode (UMI) (8bp) // [1 T] //  Insert

# After moving the UMI sequences, the script looks for a single T at the start of the sequence, and removes it.

# Script last modified 16 August 2021, Felix Krueger

polyT = {}         # storing the number of Poly Ts at the start of the read (after the UMI)
fhs = {}           # storing the filehandles for all output files

def submain():
	
	print (f"Python version: {sys.version}.")
	# allfiles = glob("*.fastq.gz")
	allfiles = sys.argv[1:]
	# print (allfiles)
	allfiles.sort() # required as glob doesn't necessarily store files in alphabetical order
	# print (allfiles)

	for filename in allfiles:
		print (f"Reading in FastQ file:\t >> {filename} <<\n")
		main(filename)
		polyT = {} # resetting

def main(filename):

	count = 0              # total sequence count
	print (f"Reading file: >{filename}<")
	
	outfh = make_out_filehandle("UMIed_plusT",filename)
	
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
			readID = readID.replace(" ","_") # this is required for e.g. Bowtie2 to retain the last part of the read ID (= the UMI sequence)

			#print (f"{readID}\n{rest}\n{line3}\n{qual_rest}\n")
			#sleep(1)

			## STEP 2: Now we remove the first T after the barcode
			# removing 1 T bp and its quality scores
			new_rest      = rest[1::]
			new_rest_qual = qual_rest[1::]
				
			# print (f"Rest:\n{rest}\n {new_rest}")
			# sleep(1)
			
			# print ("\n".join([readID, new_rest, line3, new_rest_qual]))
			outfh.write (("\n".join([readID, new_rest, line3, new_rest_qual]) + "\n").encode())
	
	outfh.close()
	
	print (f"Total number of reads processed: {count}")
	
def make_out_filehandle(sample_name,filename):
	
	print (f"Got following sample name: {sample_name}")
	
	# extracting useful parts from filename
	# Example name: lane7906_CATTTT_NCI-H23_UMI_circle_L001_R1.fastq.gz

	pattern = '(lane.*_L00\d)_(R\d.fastq.gz)'
	p = re.compile(pattern)
	print (filename)
	m = p.findall(filename)
	sample = m[0][0]
	ending = m[0][1]
	new_filename = f"{sample}_{sample_name}_{ending}"
	# print (new_filename)
	
	outfh  = gzip.open (new_filename,mode='wb',compresslevel=3)
	
	return outfh
	
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)	


if __name__ == "__main__":
	submain()
else:
	print ("Just getting imported")
