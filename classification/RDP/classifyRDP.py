#!/usr/bin/env python
# FILE: trainingresults2report.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import math
import numpy as np
import os
from Bio import SeqIO
#from keras.utils import np_utils



trainfolder = sys.argv[1]
rdpclassifierpath = sys.argv[2] # the path to the rdp classifier file classifier.jar
testfastafilename=sys.argv[3]
classificationfilename=""
if len(sys.argv) > 4:
	classificationfilename=sys.argv[4] #
classificationposition=0
if len(sys.argv) >5: 
	classificationposition=sys.argv[5] 

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(seqids,classificationfilename,pos):
	kingdompos=0
	phylumpos=1 #the position of the phyla in the classification file
	classpos=2
	orderpos=3
	familypos=4
	genuspos=5
	speciespos=6
	species=[""] * len(seqids)
	genera=[""] * len(seqids)
	families=[""] * len(seqids)
	orders=[""] * len(seqids)
	classes=[""] * len(seqids)
	phyla=[""] * len(seqids)
	kingdoms=[""] * len(seqids)
	classificationfile=open(classificationfilename)
	classifications=[""] * len(seqids)
	labels=[""] * len(seqids)
	level= 0 #kingdom
	sep="\t"
	for line in classificationfile:
		if ";" in line:
			sep=";"
		words=line.split(sep)
		if line.startswith("#"):			
			i=0
			for word in words:
				word=word.rstrip()
				if word=="species":
					speciespos=i
				if word=="genus":
					genuspos=i
				if word=="family":
					familypos=i
				if word=="order":
					orderpos=i
				if word=="class":
					orderclass=i
				if word=="phylum":
					phylumpos=i
				if word=="kingdom":
					kingdompos=i
			if pos==speciespos:
				level=6
			if pos==genuspos:
				level=5
			if pos==familypos:
				level=4
			if pos==orderpos:
				level=3
			if pos==classpos:
				level=2
			if pos==phylumpos:
				level=1	
			if pos==kingdompos:
				level=1			
			continue 
		seqid=words[0].replace(">","").rstrip()
		if seqid in seqids:
			index=seqids.index(seqid)
			classification="Root"
			kingdom="Fungi"
			if level >=0:
				if kingdompos >0:
					kingdom=words[kingdompos].rstrip()
				kingdoms[index]= kingdom
				classification=classification + ";" + kingdom
			if level >=1:
				phyla[index] = words[phylumpos].rstrip()
				classification=classification + ";" + words[phylumpos].rstrip()
			if level >=2:
				classes[index]= words[classpos].rstrip()
				classification=classification + ";" +  words[classpos].rstrip()
			if level >=3:
				orders[index]=words[orderpos].rstrip()
				classification=classification + ";" + words[orderpos].rstrip()
			if level>=4:
				families[index]=words[familypos].rstrip()
				classification=classification + ";" + words[familypos].rstrip()
			if level>=5:
				genera[index]=words[genuspos].rstrip()
				classification=classification + ";" + words[genuspos].rstrip()
			if level>=6:
				species[index]=words[speciespos].rstrip()
				classification=classification + ";" + words[speciespos].rstrip()
			classifications[index]=classification
			labels[index]=words[pos].rstrip()
	return species,genera,families,orders,classes,phyla,kingdoms,classifications,labels

def GetSeqIDs(seqrecords):
	seqids=[]
	for seqrecord in seqrecords:
		seqids.append(seqrecord.id)
	return seqids

def GetPredictedLabels(testseqids,rdpclassifiedfilename):
	predlabels =[""] * len(testseqids)
	scores=[0] * len(testseqids)
	ranks=[""] * len(testseqids)
	rdpfile=open(rdpclassifiedfilename)
	for line in  rdpfile:  
		words=line.split("\t")
		n=len(words)
		seqid = words[0]
		i=testseqids.index(seqid)
		scores[i]=float(words[n-1].rstrip())
		ranks[i]=words[n-2]
		predlabels[i]=words[n-3]
	rdpfile.close()
	return predlabels,scores,ranks

def SavePrediction(testseqids,testlabels,predlabels,probas,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tClassification\tPrediction\tProbability\n")
	for i in range(0,len(testseqids)):
		output.write(str(i) + "\t" + testseqids[i] + "\t" + testlabels[i] + "\t"  + predlabels[i] + "\t" + str(scores[i]) + "\n")
	output.close()
##############################################################################
# MAIN
##############################################################################
basefilename=GetBase(testfastafilename)
rdpclassifiedfilename = basefilename + ".rdp.out"
rdp_output=basefilename + ".rdp.classified"
#load seq records
seqrecords = list(SeqIO.parse(testfastafilename, "fasta"))
testseqids = GetSeqIDs(seqrecords)

#load given taxonomic classifcation of the test dataset
testspecies,testgenera,testfamilies,testorders,testclasses,testphyla,testkingdoms,testclassifications,testlabels=LoadClassification(testseqids,classificationfilename,int(classificationposition))

testcommand="java -Xmx1g -jar " + rdpclassifierpath + "/classifier.jar classify -t " + trainfolder + "/rRNAClassifier.properties -o " + rdpclassifiedfilename + " " + testfastafilename
os.system(testcommand)
testseqids = GetSeqIDs(seqrecords)
predlabels,scores,ranks=GetPredictedLabels(testseqids,rdpclassifiedfilename)	
SavePrediction(testseqids,testlabels,predlabels,scores,rdp_output)



