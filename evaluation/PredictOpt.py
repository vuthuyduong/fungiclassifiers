#!/usr/bin/env python
# FILE: PredictOpt.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import os
from Bio import SeqIO
import json


fastafilename = sys.argv[1] # the fasta file
classificationfilename = sys.argv[2] # the file containing taxonomic classification: class order family genus species
classificationpos = int(sys.argv[3]) # the position of the feature in the header of the sequences used for comparisons
threshold = 0.97
if "," in sys.argv[4]:
		ts = sys.argv[4].split(",")
		t=ts[len(ts)-1]
		pre_t=sys.argv[4][:len(t)]
		threshold=float(t)
else:
	threshold = float(sys.argv[4])
endthreshold = 0
if len(sys.argv) > 5:
	endthreshold = float(sys.argv[5])
step = 0.01
if len(sys.argv) > 6:
	step = float(sys.argv[6])

def GetClasses(fastafilename,classificationfilename,pos):
	seqrecords = list(SeqIO.parse(fastafilename, "fasta"))
	seqids=[]
	for seqrecord in seqrecords:
		seqids.append(seqrecord.id.split("|")[0])
	classes = []
	classnames = []
	#read the fasta file
	classificationlist= list(open(classificationfilename,"r"))
	for line  in classificationlist:
		if line.startswith("#"):
			continue
		if "," in line:
			words = line.split(',')
		else:
			words = line.split('\t')
		seqid = words[0].strip('>')
		if not seqid in seqids:
			continue
		classname = words[pos].rstrip()
		if classname in classnames:
			i=classnames.index(classname)
			refclass = classes[i] 
			refclass.append(seqid)
		else:
			classnames.append(classname)
			refclass=[]
			refclass.append(seqid)
			classes.append(refclass)
	return classes

def GetClustersByFMLC(resultfilename):
	clusters = []
	#read the clustering result file
	cluster = []
	i=0
	resultfile = open(resultfilename)
	next(resultfile)
	clusterid =-1
	for line in resultfile:	
		words = line.rstrip().split('\t')
		currentclusterid = int(words[0])
		if clusterid != currentclusterid:
			clusters.append(cluster)
			clusterid = currentclusterid 
			cluster =[]
		seqid = words[2].split('|')[0]
		cluster.append(seqid)
	clusters.append(cluster)
	return clusters


def ComputeFmeasure(classes,clusters):
	#compute F-measure
	f=0
	n=0
	for group in classes:
		m = 0
		for cluster in clusters:
			i = len(set(group) & set(cluster))
			v = float(2*i)/float((len(group) + len(cluster)))
			if m < v:
				m=v		
		n = n + len(group)
		f = f +	(len(group)*m)	
	return float(f)/float(n) 

#MAIN:
prefix = fastafilename[:-(len(fastafilename)-fastafilename.rindex("."))]
output=prefix + ".fmlc.out"
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]

if endthreshold < threshold:	
	command =  path + "./fmlc -input " + fastafilename + " -algo CCBC -thresholds " +  str(threshold) + " -output " + output
	os.system(command)
	clusters = GetClustersByFMLC(output)
	classes = GetClasses(fastafilename,classificationfilename,classificationpos)
	f = ComputeFmeasure(classes,clusters)
	print("F-measure: " + str(f))
else:
	classes = GetClasses(fastafilename,classificationfilename,classificationpos)	
	#predict optimal threshold
	t = threshold
	optthreshold=0
	maxf = 0
	outputfilename = fastafilename[:-(len(fastafilename)-fastafilename.rindex("."))] + ".opt.fmlc.out"
	resultfile = open(prefix + ".opt.fmlc.out","w")
	while t <= endthreshold:
		command =  path + "./fmlc -input " + fastafilename + " -algo CCBC -thresholds " +  str(t) + " -output fmlc_ccbc.out"
		print(command)
		os.system(command)
		clusters = GetClustersByFMLC("fmlc_ccbc.out")
		f = ComputeFmeasure(classes,clusters)
		if maxf < f:
			maxf=f
			optthreshold=t
		resultfile.write(str(t) + '\t' + str(f) + '\n')
		t = t + step
	
	resultfile.close()
	
		
