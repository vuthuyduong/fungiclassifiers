#!/usr/bin/env python
# FILE: evaluateBLAST.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import math
import numpy as np
import os
from Bio import SeqIO
#from keras.utils import np_utils
from datetime import datetime, date
from datetime import timedelta

fastafilename=sys.argv[1]
classificationfilename=sys.argv[2] #taxonomy file 
classificationlevel=int(sys.argv[3]) #the level of classification to get taxa from the taxonomy file
datapath=sys.argv[4] #training and testing datasets'path
t1=0 #the optimal threshold for classification by clustering in case t2 is not given. Otherwise it is the begin threshold to predict an optimal threshold for clustering
if len(sys.argv) >5: 
	if "," in sys.argv[5]:
		ts = sys.argv[5].split(",")
		t=ts[len(ts)-1]
		pre_t=sys.argv[5][:len(t)]
		t1=float(t)
	else:
		t1 = float(sys.argv[5])
	opthreshold=t1
t2=0 #the end threshold used to predict optimal threshold for classification
if len(sys.argv) > 6:
	t2 = float(sys.argv[6])
step = 0.001
if len(sys.argv) > 7:
	step = float(sys.argv[7])

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(classificationfilename,classificationpos):
	taxafile= open(classificationfilename)
	classification=[]
	seqnames=[]
	for line in taxafile:
		if line.startswith("#"):
			continue
		texts=line.split("\t")
		taxonname=texts[classificationlevel].rstrip()
		classification.append(taxonname)
		seqnames.append(texts[0].replace(">",""))
	return classification,seqnames

def GenerateFastafile(prefixfilename,seqids,seqrecords):	
	seqfile = open(prefixfilename + ".fasta","w")
	i=0
	for seqrecord in seqrecords:
		if i in seqids:
			#save this sequence in the seqfile
			seqfile.write(">" + seqrecord.description + "\n")
			seqfile.write(str(seqrecord.seq) + "\n")
		i=i+1	
	seqfile.close()

def PredictOpt(traindataset,classificationfilename,classificationlevel):
	if t1>=t2:
		return t1
	command=path + "PredictOpt.py " + traindataset + " " + classificationfilename + " " + str(classificationlevel) + " " + str(t1) + " " + str(t2) + " " + str(step)
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
	coverage = max(abs(maxposition - minposition) + 1,abs(refmaxposition -refminposition) + 1) + min(minposition,refminposition)-1 
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

def ComputeAccuracy(train_labels,testids,test_labels,pred_labels):
	acc=0
	identified=0
	unidentifiedlist=[]
#	acc = sum((abs(test_labels-pred_labels)>0)==False)/len(test_labels)
	for i in range(0,len(testids)):
		if test_labels[i] in train_labels:
			if test_labels[i]==pred_labels[i]:
				acc=acc+1
			identified=identified+1
		else:
			unidentifiedlist.append(testids[i])
		
	if identified >0:
		acc=float(acc)/float(identified)
	return acc,unidentifiedlist

def GetBestMatchLabels(labels,seqnames,bestmatchlist,bestscorelist,opthreshold):
	best_labels=[]
	i=0
	for seqname in bestmatchlist:
		if  seqname in seqnames:
			seqid=seqnames.index(seqname)
			if bestscorelist[i] >= optthreshold:
				best_labels.append(labels[seqid])
			else:
			#no identification
				best_labels.append(-1)
		else: #no identification
			best_labels.append(-1)
		i=i+1
	return best_labels

def GetAllLabels(trainids,train_labels,testids,test_labels):
	labels=[-1]*(len(testids)+len(trainids))
	i=0
	for testid in testids:
		labels[testid]=test_labels[i]
		i=i+1
	i=0
	for trainid in trainids:
		labels[trainid]=train_labels[i]
		i=i+1

	return labels

def SavePrediction(classification,seqnames,alllabels,testids,pred_labels,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tClassification\tPrediction\tn")
	k=0
	for i in testids:
		predictedname=""
		if pred_labels[k] in alllabels:
			j=alllabels.index(pred_labels[k])
			predictedname=classification[j]
		output.write(str(i) + "\t" + str(seqnames[i]) + "\t" + classification[i] + "\t"  + predictedname + "\n")
		k=k+1
	output.close()
##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]

reportfilename=GetBase(fastafilename) + ".blast.report"

#load seq records
seqrecords = list(SeqIO.parse(fastafilename, "fasta"))

#Load classes, classification:
classification,seqnames= LoadClassification(classificationfilename,classificationlevel)


#generate fasta files for training and test dataset
filenames = os.listdir(datapath)
prefix = GetBase(fastafilename) 
train_prefix = "train_" + prefix
test_prefix="test_"+prefix
pred_predix="preds_"+prefix

datasets=[]

for filename in filenames:
	if not ".npy" in filename:
		continue 
	basefilename = GetBase(filename)
	if (filename.startswith(test_prefix)) and ("labels" not in filename):
		testids = np.load(datapath+"/"+filename)
		GenerateFastafile(datapath+"/"+basefilename,testids,seqrecords)
	if (filename.startswith(train_prefix)) and ("labels" not in filename):
		trainids = np.load(datapath+"/"+filename)
		GenerateFastafile(datapath+"/"+basefilename,trainids,seqrecords)		
	if (filename.startswith(test_prefix)) and ("labels" in filename):
		if not basefilename.replace("test_","") in datasets:
			datasets.append(basefilename.replace("test_","").replace(".labels",""))			


#generate report
report=open(reportfilename,"w")
report.write("Test dataset\tNumber of sequences in the test dataset\tTrain dataset\tNumber of sequences in the train dataset\tOptimal threshold\tAccuracy of BLAST\tNumber of unidentified sequences by BLAST\tTraining time\tClassifying time\n")
for dataset in datasets:
	print(dataset)
	testdataset=datapath+"/test_"+dataset
	traindataset=datapath+"/train_"+dataset
	#load testids, trainids:
	testids=np.load(datapath+"/test_"+dataset+".npy")
	trainids=np.load(datapath+"/train_"+dataset+".npy")

	#load labels 
	train_labels=np.load(datapath+"/train_"+dataset+".labels.npy")
	test_labels=np.load(datapath+"/test_"+dataset+".labels.npy")
	
	#get all labels
	alllabels=GetAllLabels(trainids,train_labels,testids,test_labels)

	#predict an optimal threshold to predict a name for a sequence using Blast
	testdataset=datapath+"/test_"+dataset+".fasta"
	traindataset=datapath+"/train_"+dataset+".fasta"
	beginning=datetime.now()
	optthreshold=PredictOpt(traindataset,classificationfilename,classificationlevel)
	end=datetime.now()
	time1=(end-beginning).total_seconds()
	#search for a best match of a test sequence in a train dataset
	beginning=datetime.now()
	bestmatchlist,bestscorelist=SearchForBestMatch(testdataset,traindataset)
	#Get the best label for a test sequence based on its best match
	best_labels=GetBestMatchLabels(alllabels,seqnames,bestmatchlist,bestscorelist,optthreshold)
	end=datetime.now()
	time2=(end-beginning).total_seconds()

	#Caculate the accuracy based on searching for best match
	acc_of_searching,unidentifiedlist_by_searching=ComputeAccuracy(train_labels,testids,test_labels,best_labels)
	searching_output=datapath+"/test_"+dataset+".blast.out"

	#Save the prediction
	SavePrediction(classification,seqnames,alllabels,testids,best_labels,searching_output)

	#Print the result
	report.write(testdataset + "\t" + str(len(testids)) + "\t" + traindataset + "\t" + str(len(trainids)) +"\t" + str(optthreshold)+ "\t"+ str(acc_of_searching)+ "\t" + str(len(unidentifiedlist_by_searching)) + "\t" + str(time1) + "\t" + str(time2) + "\n")
	#Save prediction by searching
	print(testdataset + "\t" + str(len(testids)) +  "\t" + traindataset + "\t" + str(len(trainids)) + "\t" + str(optthreshold)+ "\t"+ str(acc_of_searching)+ "\t" + str(len(unidentifiedlist_by_searching)) + "\t" + str(time1) + "\t" + str(time2) + "\n")
report.close()
