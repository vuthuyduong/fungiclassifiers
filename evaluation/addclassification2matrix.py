#!/usr/bin/env python
from Bio import SeqIO
import sys
import random

matrixfilename=sys.argv[1] # the matrix obtained by fasta2matrix
classificationfilename = sys.argv[2] # the file containing  taxonomic information separated by tab
classificationposition = int(sys.argv[3]) #the position of the interest taxonomic group
output = matrixfilename[:-(len(matrixfilename)-matrixfilename.rindex("."))] + "_" + str(classificationposition) + ".matrix"
if len(sys.argv) >4:
	output = sys.argv[4]

outputFile = open(output, "w")
#read classification
classification=[]
classificationfile= open(classificationfilename)
seqids=[]
for line in classificationfile:
	if line.startswith("#"):
		continue 
	taxa=line.split("\t")
	seqids.append(taxa[0].replace(">","").rstrip())
	classification.append(taxa[classificationposition].rstrip())
matrixfile = open(matrixfilename)
header = matrixfile.readline()
outputFile.write(header)
for line in matrixfile:
	line=line.rstrip()
	seqid =line.split(",")[0]
	index=seqids.index(seqid)
	classname=classification[index]
	if classname != "":
		outputFile.write(line + ","+ classname +"\n")
outputFile.close()
