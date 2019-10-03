#!/usr/bin/env python
# FILE: classifyBLAST.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import math
import numpy as np
import os
from Bio import SeqIO
#from keras.utils import np_utils


testdataset=sys.argv[1]
traindataset=sys.argv[2]
classificationfilename = sys.argv[3]
classificationposition =int(sys.argv[4])
optthreshold=0 #the optimal threshold for classification 
if "," in sys.argv[5]:
	ts = sys.argv[5].split(",")
	t=ts[len(ts)-1]
	pre_t=sys.argv[5][:len(t)]
	optthreshold=float(t)
else:
	optthreshold = float(sys.argv[5])

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(seqIDs,classificationfilename,pos):
	classification=[""]*len(seqIDs)
	if classificationfilename == "":
		return classification
	records= list(open(classificationfilename, "r"))
	for line in records:
		if line.startswith("#"):
			continue 
		elements=line.split("\t")
		seqid = elements[0].replace(">","").rstrip()
		if seqid in seqIDs:
			index=seqIDs.index(seqid)
			classname=elements[pos].rstrip()
			classification[index]=classname
	return classification

def GetSeqIndex(seqname,seqrecords):
	i=0
	for seqrecord in seqrecords:
		if (seqname == seqrecord.id):
			return i
		i = i + 1
	return -1

def IndexSequences(filename):
	indexedfilename = GetBase(filename) + ".indexed.fasta"
	fastafile = open(filename)
	indexedfile = open(indexedfilename, "w")
	i=0
	for line in fastafile:
		if line.startswith('>'):
			indexedfile.write(">" + str(i) + "|" + line.rstrip()[1:] + "\n")
			i=i+1
		else:
			indexedfile.write(line)    
	fastafile.close()
	indexedfile.close()
	return indexedfilename

def SearchForBestMatch(testdataset,traindataset):
	#index the query sequences to be faster
	indexedtestdataset= IndexSequences(testdataset)
	#load sequeces from the fasta files
	queryrecords = list(SeqIO.parse(indexedtestdataset, "fasta"))
	bestscorelist =[0] * len(queryrecords)
	bestrefnamelist = [""] * len(queryrecords)

	#blast
	makedbcommand = "makeblastdb -in " + traindataset + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + indexedtestdataset + " -db  db -outfmt 6 -out out.txt"
	os.system(blastcommand)
	#read blast output
	blastoutputfile = open("out.txt")
	refname = ""
	queryname=""
	minposition = 0
	maxposition = 0
	refminposition = 0
	refmaxposition = 0
	identity =0
	seqlen=0
	refseqlen=0
	positions  = []
	identities =[]
	for line in blastoutputfile:
		words = line.split("\t")
		currentqueryname = words[0]
		currentrefname = words[1]
		if (refname != currentrefname) or (queryname !=currentqueryname):
			if (refname !="") and (queryname != ""):
				#compute the similarity of the previous seqquence				
				#i = GetSeqIndex(queryname, queryrecords)
				i = int(queryname.split("|")[0])
				#coverage = max(abs(maxposition - minposition) + 1,abs(refmaxposition -refminposition) + 1) + min(minposition,refminposition)-1 + min(abs(seqlen-maxposition),abs(reflen-refmaxposition))
				coverage = max(abs(maxposition - minposition) + 1,abs(refmaxposition -refminposition) + 1) + min(minposition,refminposition)-1 
				#compute identity
				identity =0 
				for pos in range(minposition,maxposition + 1):
					if pos in positions:
						indexofpos = positions.index(pos)
						identitiesatpos=identities[indexofpos]
						identity = float(identity + max(identitiesatpos)/100)
				sim=0
				if coverage > 0:
					sim  = float(identity / coverage)
				score = sim
				if  coverage < 150:
					score = float(score * coverage /  150)
				if bestscorelist[i] < score:
					bestscorelist[i]= score
					bestrefnamelist[i]=refname
			refname = currentrefname
			queryname=currentqueryname
			minposition = 0
			maxposition =0
			refminposition = 0
			refmaxposition =0 
			identity = 0
			positions = []
			identities =[]
		refpos1 =  int(words[8])
		refpos2 = int(words[9])
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		for pos in range(pos1,pos2 + 1):
			if pos in positions:
				indexofpos = positions.index(pos)
				identities[indexofpos].append(iden)
			else:
				positions.append(pos)
				identitiesatpos =[]
				identitiesatpos.append(iden)
				identities.append(identitiesatpos)
		if minposition > 0:
			minposition = min(minposition, min(pos1,pos2))
		else:
			minposition = min(pos1,pos2)
		maxposition = max(maxposition, max(pos1,pos2))
		if refminposition > 0:
			refminposition = min(refminposition, min(refpos1,refpos2))
		else:
			refminposition = min(refpos1,refpos2)
		refmaxposition = max(refmaxposition, max(refpos1,refpos2))
	#i = GetSeqIndex(queryname, queryrecords)
	i = int(queryname.split("|")[0])
	#coverage = max(maxposition - minposition + 1,refmaxposition -refminposition + 1) + min(minposition,refminposition)-1 + min(abs(seqlen-maxposition),abs(reflen-refmaxposition))
	coverage = max(maxposition - minposition + 1,refmaxposition -refminposition + 1) + min(minposition,refminposition)-1 
	#compute identity
	identity =0 
	for pos in range(minposition,maxposition + 1):
		if pos in positions:
			indexofpos = positions.index(pos)
			identitiesatpos=identities[indexofpos]
			identity = float(identity + max(identitiesatpos))/float(100)
	if coverage > 0:
		sim  = float(identity /coverage)
		score = sim
		if  coverage < 150:
			score = float(score * coverage / 150)
		if bestscorelist[i] < score:
			bestscorelist[i]= score
			bestrefnamelist[i]=refname
	os.system("rm " + indexedtestdataset)
	os.system("rm out.txt")
	return bestrefnamelist,bestscorelist

def GetBestMatchLabels(trainclassification,trainseqIDs,bestmatchlist,bestscorelist,opthreshold):
	bestlabels=[]
	i=0
	for seqid in bestmatchlist:
		if  seqid in trainseqIDs:
			index=trainseqIDs.index(seqid)
			if bestscorelist[i] >= optthreshold:
				bestlabels.append(trainclassification[index])
			else:
			#no identification
				bestlabels.append("")
		else: #no identification
			bestlabels.append("")
		i=i+1
	return bestlabels

def SavePrediction(trainclassification,testclassification,testseqIDs,pred_labels,bestscorelist,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tClassification\tPrediction\tScore\n")
	for i in range(0,len(testseqIDs)):
		testlabel=""
		if len(testclassification) > i:
			testlabel=testclassification[i]			
		output.write(str(i) + "\t" + str(testseqIDs[i]) + "\t" + testlabel + "\t"  + pred_labels[i] + "\t" + str(bestscorelist[i]) + "\n")
	output.close()

##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]
#load train seq records
trainseqrecords = list(SeqIO.parse(traindataset, "fasta"))
trainseqIDs=[]
for seq in trainseqrecords:
	trainseqIDs.append(seq.id)
#load test seq records
testseqrecords = list(SeqIO.parse(testdataset, "fasta"))
testseqIDs=[]
for seq in testseqrecords:
	testseqIDs.append(seq.id)


#Load classes, classification:
trainclassification= LoadClassification(trainseqIDs, classificationfilename, classificationposition)
testclassification= LoadClassification(testseqIDs,classificationfilename, classificationposition)
		
#search for a best match of a test sequence in a train dataset
bestmatchlist,bestscorelist=SearchForBestMatch(testdataset,traindataset)
	
#Get the best label for a test sequence based on its best match
bestlabels=GetBestMatchLabels(trainclassification,trainseqIDs,bestmatchlist,bestscorelist,optthreshold)

#Save prediction by searching 
reportfilename=GetBase(testdataset) + ".blast.classified"
SavePrediction(trainclassification,testclassification,testseqIDs,bestlabels,bestscorelist,reportfilename)
	

