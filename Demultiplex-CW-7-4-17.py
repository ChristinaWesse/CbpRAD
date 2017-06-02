#!/usr/bin/python
# Demultiplex-CW Python Script. Author: Christina Wesse (7.4.2017)
# USAGE:
# call python script with 2 arguments: 1: path/to/raw.fastq.gz 2: path/to/barcodesfile.txt
# Barcode input file must contain barcode sequences in column 1 and sample names in column 2 separated by a space
# Default restriction enzyme is KpnI ("REsite", line 73)

from __future__ import division
import gzip, sys, time, os.path											# gzip for unzipping fastq.gz, sys for user argument input, time for duration of execution, os.path for file name

def isAscii(s):															# returns True if not proper quality score ascii
	return all(ord(c) < 127 for c in s) and all(ord(c) > 32 for c in s)	# fastQ quality score start using Ascii symbols from 33 (+33 Fastq Sanger format) up to 126
		
def isDna (dna):														# returns 1 if proper DNA content
	no_bases = dna.count("A") + dna.count("G") + dna.count("C") + dna.count("T") + dna.count("N")
	if dna.count("\n")==1:												# if input str "dna" has a newline character at the end...
		if ((len(dna))-1) == no_bases:									# ... compare number of A,G,C,T and N with total length of str minus 1
			return 1
	else:																# if input str "dna" has no newline character at the end...
		if no_bases == (len(dna)):										# ... just compare number of bases with total length of str
			return 1

def checkEntry (entry):													# returns 4 if all 4 lines of one read are formated correctly
	format_ok = 0
	if entry[0].startswith("@"):										# check format of line 1. Should start with @ symbol
		format_ok = format_ok + 1																		
	if isDna(entry[1]) == 1:											# check format of line 2. Should contain DNA bases and nothing else
		format_ok = format_ok + 1
	if entry[2]==("+\n"):												# check format of line 3. Should contain a "+" symbol and a newline and nothing else
		format_ok = format_ok + 1
	if isAscii(entry[3].rstrip("\n")) is True:												# check format of line 4. Should contain ascii symbold referring to a certain quality score and nothing else
		format_ok = format_ok + 1
	return format_ok													# returns "4" if every line check was successful

def getBarcodedict (barcodefile_input):									# returns a dictionary made of the input barcode and sample names file
	dictionary = {}																
	with open(barcodefile_input,'r') as f:
		for line in f:													
			k, v  = line.split()										# separate each line by barcode ("k") and and sample name ("v")
			dictionary[k]=v												# make a dict with keys from k (barcodes) and values (sample names)
	return dictionary

singleRead = []															# global variable	
	
def readFourLines ():													# returns True as long as end of raw.fastq is not reached yet. Saves one read in global varable "singleRead"
	global singleRead
	singleRead = []	
	for i in range(4):													# read in only one single read, which contains of 4 lines within the fastq dataset
		line = f.readline()	
		if not line: 
			return False												# return False if end of raw.fastq is reached to cancel main
		else: 
			singleRead.append(line)										# make a list of one read if end of raw.fastq is not reached yet
	return True
 	
readDictionary = {}														# dictionary containing number of reads per barcode
readCount=0																# counter for total number of reads
	
def logReads ():														
	logfileName=os.path.basename(barcodefile)+".log"					# creates log file and names it after barcodesfile for identification
	f = open(logfileName, 'w')											# open file for writing
	for key in readDictionary:	
		f.write(key+":"+str(readDictionary[key])+"\n")					# write total number of reads per barcode from dictionary to log file
	f.write( "Total reads:"+str(readCount)+"\n")						# log total reads
	duration = time.time() - ticksStart									
	f.write( "Duration:" + str(duration)+"\n")							# log duration		
	f.close()

######### MAIN ##########

filename = sys.argv[1]													# raw fastq.gz file, given user input argument 1
barcodefile = sys.argv[2]												# txt file with barcode sequences and sample names, given user input argument 2
REsite = "GTACC"														# RE site pattern for KPNI
barcodesdict = getBarcodedict(barcodefile)								# make dictionary from input barcode file ("barcodefile")
																							
with gzip.open(filename, 'rb') as f:									# open temporarily the source fastq file
	ticksStart = time.time()
	while readFourLines():                      	                    # as long as end of raw.fastq.gz us not reached yet, make a list of the 4 lines (== whole read)
		readCount=readCount+1														
		if checkEntry(singleRead)!=4:									# check if the read is formatted correctly
			logfile = open("Log-Barcodes.txt", "a")						
			logfile.write("Input file is corrupted")					# log error message if read is not formatted correctly	
			continue
		dna = singleRead[1]												# the second object in the "singleRead" list is the actual dna sequence
		
		if dna.find(REsite,4)<0:										# look for cut site of restriction enzyme, but ignore first 5 bases (minimum length of barcodes; in case same sequence is within the barcode sequence)
			outputFilename="trash.fastq.gz"
			with gzip.open(outputFilename, "a") as myfile:		
				read = str(singleRead[0]) + str(singleRead[1]) + str(singleRead[2]) + str(singleRead[3]) # take the read as it is
				myfile.write(str(read))   								# save reads with no cut site in this file ("trash")			
				readDictionary[outputFilename] = readDictionary.get(outputFilename,0) +1	# set counter + 1 for "trash" reads			
				continue
														
		start = dna.find(REsite,4)										# store start position of cut site			
		brcd = dna[:start]												# set everything in "front" of the cut site as barcode (brcd) 
		if brcd not in barcodesdict:									# if "brcd" does not match exactly with a key within the barcode dict log error message 		
			outputFilename="no_barcodes.fastq.gz"
			with gzip.open(outputFilename, "a") as myfile:				# append to log file
				read = str(singleRead[0]) + str(singleRead[1]) + str(singleRead[2]) + str(singleRead[3]) 	# save 4 lines of read to variable
				myfile.write(str(read))									# save variable with 4 lines to logfile	
				readDictionary[outputFilename] = readDictionary.get(outputFilename,0) +1  #increase counter for no_barcode found
				continue
					
		samplename = barcodesdict.get(brcd)								# ... get the according sample name from dictionay	
		outputFilename="sample_" + str(samplename) + ".fastq.gz"		# create outputfilename with name of the sample
		with gzip.open(outputFilename, "a") as myfile:	
			read = str(singleRead[0]) + str((singleRead[1])[(len(brcd)):]) + str(singleRead[2]) + str((singleRead[3])[(len(brcd)):])	# trim barcode sequence in line two and quality score score in line four during this step 
			myfile.write(str(read))										# save read with proper cut site and barcode in separate fastq.gz file	
			readDictionary[outputFilename] = readDictionary.get(outputFilename,0) +1   #increase counter for sample
		if readCount % 100000 == 0:										# only log every 100.000 read
			logReads()			
	logReads()															# Log at end of program
