#!/usr/bin/env python
import sys, gzip, os
import getopt
from glob import glob
from time import sleep
import re

# This script removes the first 8bp of a read and adds the sequence to the readID. The quality information is discarded
# In addition, the sample level barcode will be written into the filename (quality information is also discarded)

# This is the expected structure of the FastQ files will be:

# NNNNNNNNBBBB(A)nSEQUENCESPECIFIC whereas before it was NNNNNNNN(A)nSEQUENCESPECIFIC
# NNNNNNNN (UMI)(8bp) // 
# BBBB is the sample barcode (currently either AGTC or GACT) //
# PolyA (A)n is the poly(A) //      
# Insert

# After moving the UMI and sample-level barcode sequences, the script looks for up to 3 T at the start of the sequence, and removes those.
# Sequences with more than 3 Ts at the 5' end are clipped a maximum of 3 TTT

# Script last modified  21 February 2022 to add some additional indexes, Felix Krueger

polyT = {}         # storing the number of Poly Ts at the start of the read (after the UMI)
fhs = {}           # storing the filehandles for all output files

def submain():
	
	print (f"Python version: {sys.version}.")
	# allfiles = glob("*.fastq.gz")
	# Changing to command line given list of arguments as infiles
	allfiles = sys.argv[1:]
	# print (allfiles)

	allfiles.sort() # required as glob doesn't necessarily store files in alphabetical order
	# print (allfiles)

	for filename in allfiles:
		print (f"Reading in FastQ file:\t >> {filename} <<\n")
		main(filename)
		polyT.clear() # resetting
		fhs.clear()   # resetting

def main(filename):

	expected_count   = 0
	unexpected_count = 0   
	count = 0              # total sequence count
	
	print (f"Reading file: >{filename}<")
	
	faithless_barcodes = {}

	# making output filehandles
	make_out_filehandle("UMIed",filename)
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
			# barcode (8bp)    //   sample-level barcode     //    PolyT     //      Insert
		
			barcode   = seq[0:8]
			sampleBarcode = seq[8:12:]
			rest      = seq[12::]
			qual_rest = qual[12::]
			# print (f"sequence: {seq}\numi:  {barcode}\nsample-level barcode:    {sampleBarcode}\nrest:             {rest}\n")
			readID += f':{sampleBarcode}:{barcode}'

			#print (f"{readID}\n{rest}\n{line3}\n{qual_rest}\n")
			#sleep(1)

			## STEP 2: Now we need to find a number of PolyTs at the start

			pattern = '^T+'
			p = re.compile(pattern)
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

			### First generation

			# sample_level_barcode_1 = "AGTC"  ### first generation - now with a new supplier (Sept 03, 2021)
			# # sample_level_barcode_2 = "GACT"  ### first generation
			# sample_level_barcode_3 = "CTTG" 
			# # sample_level_barcode_4 = "TCGA"
			# sample_level_barcode_5 = "AAGG"
			# sample_level_barcode_6 = "TTCC" 
			# sample_level_barcode_7 = "GTGC" 
			# sample_level_barcode_8 = "GCCA" 
			# sample_level_barcode_9 = "GATG" 

			# currently indexes 1, 3, 5, 6, 7, 8, 9
			if sampleBarcode == "AGTC" or sampleBarcode == "CTTG" or sampleBarcode == "TCGA" or sampleBarcode == "TTCC" or sampleBarcode == "AAGG"or sampleBarcode == "GTGC" or sampleBarcode == "GCCA" or sampleBarcode == "GATG":
				# print (f"Expected: {sampleBarcode}")
				fhs[sampleBarcode].write (("\n".join([readID, new_rest, line3, new_rest_qual]) + "\n").encode())
			else:
				# print (f"Unexpected/Unwanted")
				fhs["unassigned"].write (("\n".join([readID, new_rest, line3, new_rest_qual]) + "\n").encode())
				if sampleBarcode not in faithless_barcodes.keys():
					faithless_barcodes[sampleBarcode] = 0
				faithless_barcodes[sampleBarcode] += 1
				# pass
			# sleep(1)
			
	
	close_filehandles()
	
	barcode_count = 0
	for rogue in sorted (faithless_barcodes, key=faithless_barcodes.get, reverse=True):
		print (f"{rogue}\t{faithless_barcodes[rogue]}")
		barcode_count += 1
		if barcode_count == 50:
			break

	print (f"Total number of reads processed: {count}")
	t_count = 0
	for tees in sorted (polyT.keys()):
		print (f"{tees}\t{polyT[tees]}")
		t_count += 1
		if t_count == 10:
			break

def make_out_filehandle(sample_name,filename):
	
	print (f"Got following sample name: {sample_name}")
	
	# extracting useful parts from filename
	# Example name: lane7265_ACTTGA_fob1_YPD_LIGseq_L001_R1.fastq.gz
	
	### Update 21 July 2021
	# 	TrAEL index	In the adaptor	What is actually read
	# 1	GACT	agtc - no longer used
	# 2	AGTC	gact - no longer used
	# 3	CAAG	cttg
	# 4	TCGA	tcga
	# 5	GGAA	ttcc
	# 6	CCTT	aagg

	### Update 21 February 2022
	# Index	In the adaptor	What is actually read	
	# 1	GACT	agtc	
	# 2	AGTC	gact	
	# 3	CAAG	cttg	
	# 4	TCGA	tcga	
	# 5	CCTT	aagg	
	# 6	GGAA	ttcc	new
	# 7	GCAC	gtgc	new
	# 8	TGGC	gcca	new
	# 9	CATC	gatg	new

	# We will also need to add the sample level barcodes to the filename.
	sample_level_barcode_1 = "AGTC"  ### first generation - now with a new supplier (Sept 03, 2021)
	sample_level_barcode_2 = "GACT"  ### first generation
	sample_level_barcode_3 = "CTTG" 
	sample_level_barcode_4 = "TCGA"
	sample_level_barcode_5 = "AAGG"
	sample_level_barcode_6 = "TTCC" 
	sample_level_barcode_7 = "GTGC" 
	sample_level_barcode_8 = "GCCA" 
	sample_level_barcode_9 = "GATG" 
	
	pattern = '(lane.*_L00\d)_(R\d.fastq.gz)'
	p = re.compile(pattern)
	print (filename)
	m = p.findall(filename)
	sample = m[0][0]
	ending = m[0][1]

	new_filenames = []
	
	new_filename_1 = f"{sample}_{sample_name}_{sample_level_barcode_1}_index1_{ending}"
	new_filenames.append(f"{new_filename_1}:{sample_level_barcode_1}")
	
	new_filename_2 = f"{sample}_{sample_name}_{sample_level_barcode_2}_index2_{ending}"
	new_filenames.append(f"{new_filename_2}:{sample_level_barcode_2}")
	
	new_filename_3 = f"{sample}_{sample_name}_{sample_level_barcode_3}_index3_{ending}"
	new_filenames.append(f"{new_filename_3}:{sample_level_barcode_3}")
	
	new_filename_4 = f"{sample}_{sample_name}_{sample_level_barcode_4}_index4_{ending}"
	new_filenames.append(f"{new_filename_4}:{sample_level_barcode_4}")
	
	new_filename_5 = f"{sample}_{sample_name}_{sample_level_barcode_5}_index5_{ending}"
	new_filenames.append(f"{new_filename_5}:{sample_level_barcode_5}")
	
	new_filename_6 = f"{sample}_{sample_name}_{sample_level_barcode_6}_index6_{ending}"
	new_filenames.append(f"{new_filename_6}:{sample_level_barcode_6}")

	new_filename_7 = f"{sample}_{sample_name}_{sample_level_barcode_7}_index7_{ending}"
	new_filenames.append(f"{new_filename_7}:{sample_level_barcode_7}")

	new_filename_8 = f"{sample}_{sample_name}_{sample_level_barcode_8}_index8_{ending}"
	new_filenames.append(f"{new_filename_8}:{sample_level_barcode_8}")

	new_filename_9 = f"{sample}_{sample_name}_{sample_level_barcode_9}_index9_{ending}"
	new_filenames.append(f"{new_filename_9}:{sample_level_barcode_9}")

	# Unassigned file
	new_filename_10 = f"{sample}_{sample_name}_unassigned_{ending}"
	new_filenames.append(f"{new_filename_10}:unassigned")

	for new_fh in new_filenames:
		open_filehandles(new_fh.split(":")[1], new_fh.split(":")[0])

def open_filehandles(sample_level_barcode, fname):
		print (f"Opening filehandle for {sample_level_barcode} and {fname}")
		fhs[sample_level_barcode] = gzip.open (fname,mode='wb',compresslevel=3)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)	

def close_filehandles():
	for name in fhs.keys():
		fhs[name].close()

if __name__ == "__main__":
	submain()
else:
	print ("Just getting imported")
