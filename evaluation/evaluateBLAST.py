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
from sklearn.metrics import precision_recall_fscore_support
from datetime import datetime, date
from datetime import timedelta
import multiprocessing
import random

nproc=multiprocessing.cpu_count()

fastafilename=sys.argv[1]
classificationfilename=sys.argv[2] #taxonomy file 
classificationlevel=int(sys.argv[3]) #the level of classification to get taxa from the taxonomy file
datapath=sys.argv[4] #training and testing datasets'path
t1=0 #the optimal threshold for classification by clustering in case t2 is not given. Otherwise it is the begin threshold to predict an optimal threshold for clustering
if len(sys.argv) >5: 
	t1 = float(sys.argv[5])
opthreshold=t1
t2=1 #the end threshold used to predict optimal threshold for classification
if len(sys.argv) > 6:
	t2 = float(sys.argv[6])
step = 0.001
if len(sys.argv) > 7:
	step = float(sys.argv[7])
mincoverage=300
if len(sys.argv) > 8:
	mincoverage = int(sys.argv[8])
reportfilename=""
if len(sys.argv) > 9:
	reportfilename=sys.argv[9]
evalnumber=0
if len(sys.argv) >10:
	evalnumber=int(sys.argv[10])	

def GetBase(filename):
	base=filename
	if "." in filename:
		base=filename[:-(len(filename)-filename.rindex("."))]
	if "/" in base:
		base=base[base.rindex("/")+1:]
	return base 

class PointDef:
	def __init__(self, id,name,flag,neighbors,seqrecord):
		self.id=id
		self.name=name
		self.flag=flag
		self.neighbors = neighbors
		self.seqrecord=seqrecord

class ClusterDef:
	def __init__(self, id, pointids):
		self.id=id
		self.pointids=pointids

def LoadNeighbors(seqrecords,scorematrix,threshold):
	neighborlist=[]
	for i in range(len(seqrecords)):
		neighborlist.append([])
	for i in range(0,len(seqrecords)-1):
		for j in range(i+1,len(seqrecords)):
				if scorematrix[i][j] >= threshold:
					if (j not in neighborlist[i]):
						neighborlist[i].append(j)
					if (i not in neighborlist[j]):
						neighborlist[j].append(i)
	#os.system("rm out.txt")
	return neighborlist

def LoadPoints(neigborlist,seqrecords):
	points=[]
	i=0
	for seqrecord in  seqrecords:
		point = PointDef(i,seqrecord.id,False,neigborlist[i],seqrecord)
		points.append(point)
		i=i+1
	return points

def ExpandCluster(root,cluster,points):
	for i in root.neighbors:
		if points[i].flag==False:
			points[i].flag=True
			cluster.pointids.append(i)
			#print("cluster id:" + str(cluster.id) + "\t" + str(len(cluster.pointids)))
			ExpandCluster(points[i],cluster,points)		

def Cluster(points,clusters):
	for point in points:
		if point.flag==False:
			point.flag=True
			cluster = ClusterDef(len(clusters),[])
			cluster.pointids.append(point.id)
			#print("cluster id:" + str(cluster.id) + "\t" + str(len(cluster.pointids)))
			ExpandCluster(point, cluster,points)
			clusters.append(cluster)
		
def ComputeFmeasure(classes,clusters):
	#compute F-measure
	f=0
	n=0
	for group in classes:
		m = 0
		for cluster in clusters:
			i = len(set(group) & set(cluster.pointids))
			v = float(2*i)/float((len(group) + len(cluster.pointids)))
			if m < v:
				m=v		
		n = n + len(group)
		f = f +	(len(group)*m)	
	return float(f)/float(n) 

def LoadClasses(seqrecords,classificationfilename,pos):
	#records= open(classificationfilename,errors='ignore')
	records= open(classificationfilename)
	allseqids=[]
	allclassification=[]
	for line in records:
		words=line.split("\t")
		seqid = words[0].rstrip().replace(".1","").replace(">","")
		if pos < len(words):
			allseqids.append(seqid)
			allclassification.append(words[pos].rstrip())
	classes = []
	classnames = []
	classification=[]
	i=0
	for seqrecord in seqrecords:
		classname=""
		if seqrecord.id in allseqids:
			index=allseqids.index(seqrecord.id)
			classname=allclassification[index]
		if classname in classnames:
			classid=classnames.index(classname)
			classes[classid].append(i)
		else:
			classnames.append(classname)
			refclass=[]
			refclass.append(i)
			classes.append(refclass)
		classification.append(classname)
		i=i+1
	return classes,classnames,classification

def LoadClassification(seqrecids,classificationfilename,classificationpos):
	taxafile= open(classificationfilename)
	classification=[""]*len(seqrecids)
	for line in taxafile:
		if line.startswith("#"):
			continue
		texts=line.split("\t")
		seqid=texts[0].replace(">","")
		if not seqid in seqrecids:
			continue
		taxonname=""
		if classificationlevel < len(texts):
			taxonname=texts[classificationlevel].rstrip()
		classification[seqrecids.index(seqid)]=taxonname
	return classification

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

def PredictOpt(traindataset,classificationfilename,classificationpos,t1,t2,step,mincoverage):
	optfilename=traindataset[:-(len(traindataset)-traindataset.rindex("."))] + ".opt"
	optthreshold=t1
	bestFmeasure=0
	if os.path.isfile(optfilename):
		with open(optfilename, 'r') as f:
			lines = f.read().splitlines()  
			lastline = lines[-1]
		optthreshold=float((lastline.split("BestFmeasure:")[0]).split("Opt threshold:")[1])
		bestFmeasure=float(lastline.split("BestFmeasure:")[1])
		return optthreshold,bestFmeasure
	#load sequences
	seqrecords=list(SeqIO.parse(traindataset, "fasta"))
	#load classification at given positions 
	classes,classnames,classification=LoadClasses(seqrecords,classificationfilename,classificationpos)
	sys.setrecursionlimit(len(seqrecords)*2)
	#seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load similarity matrix
	scorematrix=ComputeBLASTscoreMatrix(traindataset,seqrecords,mincoverage)
	t=t1
	#compute optimal threshold
	optfile=open(optfilename,"w")
	optfile.write("Threshold\tFmeasure\n")
	while t < t2:
		#load neighbors
		neighborlist = LoadNeighbors(seqrecords,scorematrix,t)	
		points=LoadPoints(neighborlist,seqrecords)
		clusters=[]
		Cluster(points,clusters)	
		fmeasure=ComputeFmeasure(classes,clusters)
		optfile.write(str(t) + "\t" + str(fmeasure) + "\n" )
		if fmeasure > bestFmeasure:
			bestFmeasure=fmeasure
			optthreshold=t
		t=t+step	
	optfile.write("Opt threshold: " + str(optthreshold) + "\tBestFmeasure: " +  str(bestFmeasure) + "\n")
	optfile.close()
	return optthreshold,bestFmeasure

def ComputeVariation(classname,trainclass,trainrecords,mincoverage):
	records=[]
	for i in trainclass:
		records.append(trainrecords[i])
	#save records to fastafile
	trainclassfilename=classname.replace(" ","_") + ".fasta"
	SeqIO.write(records, trainclassfilename, "fasta")
	threshold=1
	scorematrix=ComputeBLASTscoreMatrix(trainclassfilename,records,mincoverage)
	scorelist=[]
	for i in range(0,len(scorematrix)-2):
		for j in range(i+1,len(scorematrix)-1):
			scorelist.append(scorematrix[i][j])
	threshold=1
	minthreshold=1		
	if len(scorelist) >0:
		x=np.array(scorelist)
		minthreshold=min(x)
		threshold=np.median(x)
#		if threshold < np.mean(x)-0.2:
	if len(scorelist) >2:		
		limit=np.percentile(x,50)
		filtered_x=x[x>=limit]
		threshold=np.median(filtered_x)
		#threshold=statistics.median(scorelist)
	os.system("rm " + trainclassfilename)	
	return threshold,minthreshold

def EvaluateVariation(classname,trainclass,trainrecords,mincoverage):
	thresholds=[]
	minthresholds=[] 
	for i in range(0,10):
		n=int(len(trainclass)/10)
		selectedindexes=random.sample(range(0, len(trainclass)), k=n)
		selectedsequences=[]
		for index in selectedindexes:
			selectedsequences.append(trainclass[index])
		threshold,minthreshold=ComputeVariation(classname,selectedsequences,trainrecords,mincoverage)
		thresholds.append(threshold)
		minthresholds.append(minthreshold)
	threshold = np.median(np.array(thresholds))
	minthreshold =  np.median(np.array(minthresholds))
	return threshold,minthreshold

def ComputeVariations(traindataset,classificationfilename,classificationpos,mincoverage):
	filename=traindataset[:-(len(traindataset)-traindataset.rindex("."))] + ".variation"
	if os.path.isfile(filename):
		trainclassnames,thresholds=LoadVariations(traindataset)
		return trainclassnames,thresholds
	#load train records
	trainrecords = list(SeqIO.parse(traindataset, "fasta"))
	#load classification at given positions 
	trainclasses,trainclassnames,trainclassification=LoadClasses(trainrecords,classificationfilename,classificationpos)
	i=0
	thresholds=[]
	variationfile=open(filename,"w")
	variationfile.write("Classname\tThreshold\tMin. threshold\tNumber of the sequences\n")
	for trainclass in trainclasses:
		threshold=1
		minthreshold=1
		if len(trainclass) <100:
			threshold,minthreshold=ComputeVariation(trainclassnames[i],trainclass,trainrecords,mincoverage)
		else:
			threshold,minthreshold=EvaluateVariation(trainclassnames[i],trainclass,trainrecords,mincoverage)
		thresholds.append(threshold)
		variationfile.write(trainclassnames[i] + "\t" + str(threshold) + "\t" + str(minthreshold) + "\t" + str(len(trainclass)) + "\n")
		i=i+1
	variationfile.close()
	return trainclassnames,thresholds

def LoadVariations(traindataset):
	filename=traindataset[:-(len(traindataset)-traindataset.rindex("."))] + ".variation"
	variationfile=open(filename)
	next(variationfile)
	trainclassnames=[]
	thresholds=[]
	for line in variationfile:
		words=line.rstrip().split("\t")
		classname=words[0]
		threshold=float(words[1])
		trainclassnames.append(classname)
		thresholds.append(threshold)
	variationfile.close()
	return trainclassnames,thresholds

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

def ComputeBLASTscoreMatrix(fastafilename,seqrecords,mincoverage):
	scorematrix = [[0 for x in range(len(seqrecords))] for y in range(len(seqrecords))] 
	seqids=[]
	for rec in seqrecords:
		seqids.append(rec.id)
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	refid = ""
	score=0
	queryid=""
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
		queryid = words[0].rstrip()
		refid = words[1].rstrip()
		i = seqids.index(queryid)
		j = seqids.index(refid)
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage
		if scorematrix[i][j] < score:
			scorematrix[i][j]=score
			scorematrix[j][i]=score
	os.system("rm out.txt")
	return scorematrix


def ComputeBestBLASTscore(query,reference,mincoverage):
	indexed_query= IndexSequences(query)

	#load sequeces from the fasta files
	queryrecords = list(SeqIO.parse(indexed_query, "fasta"))
	#refrecords = list(SeqIO.parse(reference, "fasta"))

	bestscorelist =[0] * len(queryrecords)
	bestrefidlist = [""] * len(queryrecords)

	#blast
	makedbcommand = "makeblastdb -in " + reference + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + indexed_query + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	#blastcommand = "blastn -query " + indexed_query + " -subject " + reference + " -outfmt 6 -out out.txt"
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	refid = ""
	
	score=0
	queryid=""
	for line in blastoutputfile:
		words = line.split("\t")
		queryid=words[0]
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		refid=words[1]
		score=sim
		if coverage < mincoverage:
				score=float(score * coverage)/mincoverage
		i = int(queryid.split("|")[0])		
		if score > bestscorelist[i]:
			bestscorelist[i]= score
			bestrefidlist[i]=refid
	os.system("rm " + indexed_query)		
	os.system("rm out.txt")
	return bestrefidlist,bestscorelist
	
def CalculateMetrics(train_labels,testids,test_labels,pred_labels): 
	precision,recall,fscore,support=precision_recall_fscore_support(test_labels,pred_labels,average='macro')
	return precision,recall,fscore

def GetBestMatchLabels(bestmatchlist,bestscorelist,labels,seqrecids,classification,classnames,thresholds,optthreshold):
	best_labels=[]
	i=0
	for seqrecid in bestmatchlist:
		if  seqrecid in seqrecids:
			opt=optthreshold
			if opt < 0:
				classname=classification[seqrecids.index(seqrecid)]
				threshold=thresholds[classnames.index(classname)]
				opt=threshold
			if bestscorelist[i] >= opt:
				best_labels.append(labels[seqrecids.index(seqrecid)])
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

def SavePrediction(classification,seqrecids,labels,testids,pred_labels,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tClassification\tPrediction\n")
	i=0
	for testid in testids:
		predictedname=""
		if pred_labels[i] in labels:
			j=labels.index(pred_labels[i])
			predictedname=classification[j]
		output.write(str(testid) + "\t" + str(seqrecids[testid]) + "\t" + classification[testid] + "\t"  + predictedname + "\n")
		i=i+1
	output.close()
##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]

if reportfilename =="":
	reportfilename=GetBase(fastafilename) + ".blast.report"

#load seq records
seqrecords = list(SeqIO.parse(fastafilename, "fasta"))
seqrecids=[]
for seqrec in seqrecords:
	seqrecids.append(seqrec.id)

#Load classes, classification:
classification= LoadClassification(seqrecids,classificationfilename,classificationlevel)
#generate fasta files for training and test dataset
filenames = os.listdir(datapath)
prefix = GetBase(fastafilename) 
if "." in prefix:
	prefix=prefix[0:prefix.index(".")]
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

#check if the output file exists
reportexists=False
evaluateddatasets=[]
if os.path.isfile(reportfilename)==True:
	report=open(reportfilename)
	next(report)
	for line in report:
		if line.rstrip()=="":
			continue
		testdataset=line.split("\t")[0]
		dataset=GetBase(testdataset).split("test_")[1]
		print(dataset)
		evaluateddatasets.append(dataset)
	report.close()
	reportexists=True

#generate report
if reportexists==False:
	report=open(reportfilename,"w")
	report.write("Test dataset\tNumber of sequences in the test dataset\tTrain dataset\tNumber of sequences in the train dataset\tNumber of unidentified sequences\tTraining time\tClassifying time\tOptimal threshold\tBest Fmeasure\n")
	report.close()
i=0
for dataset in datasets:
	if dataset in evaluateddatasets:
		continue
	print(dataset)
#	testdataset=datapath+"/test_"+dataset
#	traindataset=datapath+"/train_"+dataset
	#load testids, trainids:
	testids=np.load(datapath+"/test_"+dataset+".npy")
	trainids=np.load(datapath+"/train_"+dataset+".npy")

	#load labels 
	train_labels=np.load(datapath+"/train_"+dataset+".labels.npy")
	test_labels=np.load(datapath+"/test_"+dataset+".labels.npy")
	
	#get all labels
	labels=GetAllLabels(trainids,train_labels,testids,test_labels)

	#predict an optimal threshold to predict a name for a sequence using Blast
	testdataset=datapath+"/test_"+dataset+".fasta"
	traindataset=datapath+"/train_"+dataset+".fasta"
	beginning=datetime.now()
	#optthreshold=PredictOpt(traindataset,classificationfilename,classificationlevel)
	optthreshold,bestFmeasure=PredictOpt(traindataset,classificationfilename,classificationlevel,t1,t2,step,mincoverage)
	#compute median threshold for each of the group of the training dataset
	trainclassnames=[]
	thresholds=[]
	if optthreshold==0:
		trainclassnames,thresholds=ComputeVariations(traindataset,classificationfilename,classificationlevel,mincoverage)
	end=datetime.now()
	time1=(end-beginning).total_seconds()
	#search for a best match of a test sequence in a train dataset
	beginning=datetime.now()
	#bestmatchlist,bestscorelist=SearchForBestMatch(testdataset,traindataset)
	bestmatchlist,bestscorelist=ComputeBestBLASTscore(testdataset,traindataset,mincoverage)
	#Get the best label for a test sequence based on its best match
	best_labels=GetBestMatchLabels(bestmatchlist,bestscorelist,labels,seqrecids,classification,trainclassnames,thresholds,optthreshold)
	end=datetime.now()
	time2=(end-beginning).total_seconds()
	
	#calulate metrics
	#precision,recall,fscore,unidentifiedlist=CalculateMetrics(train_labels,testids,test_labels,best_labels)
	
	#get unidentified list
	unidentifiedlist=[]
	for i in range(0,len(testids)):
		if test_labels[i] not in train_labels:
			unidentifiedlist.append(testids[i])
			
	searching_output=datapath + "/test_" + dataset + ".blast.classified"
	
	#Save the prediction
	SavePrediction(classification,seqrecids,labels,testids,best_labels,searching_output)
	#Print the result
	report=open(reportfilename,"a")
	report.write(GetBase(testdataset) + "\t" + str(len(testids)) + "\t" + GetBase(traindataset) + "\t" + str(len(trainids))  + "\t" + str(len(unidentifiedlist)) + "\t" + str(time1) + "\t" + str(time2) + "\t" + str(optthreshold) + "\t" + str(bestFmeasure) + "\n")
	report.close()
	print(GetBase(testdataset) + "\t" + str(len(testids)) +  "\t" + GetBase(traindataset) + "\t" + str(len(trainids)) + "\t" + str(len(unidentifiedlist)) + "\t" + str(time1) + "\t" + str(time2) + "\t" + str(optthreshold) + "\t" + str(bestFmeasure) +  "\n")
	i=i+1
	if i==evalnumber:
		break
print("The running time of the model is saved in the file " + reportfilename + ". The results of the classification are saved as .blast.classified files in the folder " + datapath + "." )
