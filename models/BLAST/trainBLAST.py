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

#looking for the best threshold for classification
traindataset=sys.argv[1]
trainclassificationfilename = sys.argv[2]
trainclassificationposition =int(sys.argv[3])
t1=0 #the optimal threshold for classification by clustering in case t2 is not given. Otherwise it is the begin threshold to predict an optimal threshold for clustering 
if "," in sys.argv[4]:
	ts = sys.argv[4].split(",")
	t=ts[len(ts)-1]
	pre_t=sys.argv[4][:len(t)]
	t1=float(t)
else:
	t1 = float(sys.argv[4])
opthreshold=t1
t2=0 #the end threshold used to predict optimal threshold for classification
t2 = float(sys.argv[5])
step = float(sys.argv[6])
if step==0:
	step=0.001

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

def PredictOpt(traindataset,trainclassificationfilename,trainclassificationpos,t1,t2,step):
	if t1>=t2:
		return t1
	command=path + "PredictOpt.py " + traindataset + " " + trainclassificationfilename + " " + str(trainclassificationposition) + " " + str(t1) + " " + str(t2) + " " + str(step)
	os.system(command)
	optfilename=GetBase(traindataset) + ".opt.fmlc.out"
	optfile= open(optfilename)
	thresholds=[]
	scores=[]
	for line in optfile:
		t=float(line.rstrip().split("\t")[0])
		s=float(line.rstrip().split("\t")[1])
		scores.append(s)
		thresholds.append(t)
	bestscore=max(scores)
	index=scores.index(bestscore)
	optthreshold=thresholds[index]

	return optthreshold


##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]

#load seq records
trainseqrecords = list(SeqIO.parse(traindataset, "fasta"))
trainseqIDs=[]
for seq in trainseqrecords:
	trainseqIDs.append(seq.id)

#Load classes, classification:
trainclassification= LoadClassification(trainseqIDs, trainclassificationfilename, trainclassificationposition)
	
#predict optimal threshold	
optthreshold=t1
if t1 < t2:
	optthreshold=PredictOpt(traindataset,trainclassificationfilename,trainclassificationposition,t1,t2,step)

print("The optimal threshold for classification: " + str(optthreshold))
	
	

