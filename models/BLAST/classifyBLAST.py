#!/usr/bin/env python
# FILE: classifyBLAST.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
if sys.version_info[0] >= 3:
	unicode = str
import numpy as np
import os, argparse
from Bio import SeqIO
import json
import multiprocessing
import random

nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='classifyBLAST.py',  
							   usage="%(prog)s [options] -i fastafile -r referencefastafile -t optthreshold -cn classificationfilename -p classificationposition -mp minproba -mc mincoverage ",
							   description='''Script that classifies the sequences of the fasta files using BLAST classification. The classified sequences with a probability less than minproba will be verified using BLAST. The mincoverage is given for BLAST comparision. The json file is the the variation of the sequences within each group of the training dataset. This file is created during the training of the model, and used optionally for the verification of the classification. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be classified.')
parser.add_argument('-r','--reference', required=True, help='the reference fasta file.')
parser.add_argument('-o','--out', help='The folder name containing the model and associated files.') #optional
parser.add_argument('-cn','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-p','--classificationpos', required=True, type=int, default=0, help='the classification position to load the classification.')
parser.add_argument('-t','--threshold', required=True, type=float, default=0.97, help='The threshold for the classification.')
parser.add_argument('-mp','--minproba', type=float, default=1.1, help='Optional. The minimum probability for verifying the classification results.')
parser.add_argument('-mc','--mincoverage', type=int, default=300, help='Optinal. Minimum coverage required for the identitiy of the BLAST comparison.')
parser.add_argument('-j','--variation', help='Optinal. The json file containing the variation within each class of the classifier.')

args=parser.parse_args()
testdataset= args.input
traindataset = args.reference
optthreshold=args.threshold
minprobaforBlast=args.minproba
mincoverage = args.mincoverage
classificationfilename=args.classification
classificationposition=args.classificationpos
jsonvariationfilename=args.variation
reportfilename=args.out

#testdataset=sys.argv[1]
#traindataset=sys.argv[2]
#classificationfilename = sys.argv[3]
#classificationposition =int(sys.argv[4])
#mincoverage=300
#if len(sys.argv) >5:
#	mincoverage = float(sys.argv[5])
#optthreshold=-1 #the optimal threshold for classification, if optthreshold=-1, the threshold of the group will be taken 
#if len(sys.argv) >6:
#	optthreshold = float(sys.argv[6])	
#jsonvariationfilename = "" #json format, optional
#if len(sys.argv) >7:
#	jsonvariationfilename = sys.argv[7]
def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(seqIDs,seqrecords,classificationfilename,pos):
	classification=[""]*len(seqIDs)
	classes=[]
	classnames=[]
	level=str(pos)
	if classificationfilename == "":
		return classification
	records= list(open(classificationfilename, "r"))
	for line in records:
		elements=line.split("\t")
		if line.startswith("#"):
			level=elements[pos]
			continue 
		seqid = elements[0].replace(">","").rstrip()
		if seqid in seqIDs:
			index=seqIDs.index(seqid)
			classname=elements[pos].rstrip()
			classification[index]=classname
			if classname in classnames:
				classid=classnames.index(classname)
				classes[classid].append(seqrecords[index])
			else:
				classnames.append(classname)
				seqs=[]
				seqs.append(seqrecords[index])
				classes.append(seqs)
	return classification,classes,classnames,level

def GetSeqIndex(seqname,seqrecords):
	i=0
	for seqrecord in seqrecords:
		if (seqname == seqrecord.id):
			return i
		i = i + 1
	return -1

def ComputeBLASTscoreMatrix(fastafilename,records,mincoverage):
	scorematrix = [[0 for x in range(len(records))] for y in range(len(records))] 
	seqids=[]
	for rec in records:
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

def ComputeVariation(reffilename,mincoverage):
	#load sequeces from the fasta files
	records = list(SeqIO.parse(reffilename, "fasta"))
	scorematrix=ComputeBLASTscoreMatrix(reffilename,records,mincoverage)
	scorelist=[]
	for i in range(0,len(scorematrix)-2):
		for j in range(i+1,len(scorematrix)-1):
			if i!=j:
				scorelist.append(scorematrix[i][j])
	threshold=1
	minthreshold=1		
	if len(scorelist) >0:
		x=np.array(scorelist)
		threshold=np.median(x)
		minthreshold=np.min(x)
	return threshold,minthreshold

def EvaluateVariation(taxonname,sequences,mincoverage):
	thresholds=[]
	minthresholds=[] 
	for i in range(0,10):
		n=int(len(sequences)/10)
		selectedindexes=random.sample(range(0, len(sequences)), k=n)
		selectedsequences=[]
		for index in selectedindexes:
			selectedsequences.append(sequences[index])
		fastafilename=taxonname.replace(" ","_") + ".fasta"
		SeqIO.write(selectedsequences,fastafilename,"fasta")
		threshold,minthreshold=ComputeVariation(fastafilename,mincoverage)
		os.system("rm " + fastafilename)
		thresholds.append(threshold)
		minthresholds.append(minthreshold)
	threshold = np.median(np.array(thresholds))
	minthreshold =  np.median(np.array(minthresholds))
	return threshold,minthreshold

def ComputeVariations(variationfilename,classes,classnames,mincoverage):
	#create json dict
	variations={}
	i=0
	for taxonname in classnames:
		sequences=classes[i]
		if len(sequences) >0:
			threshold=0
			minthreshold=0
			if len(sequences) < 100:				
				fastafilename=taxonname.replace(" ","_") + ".fasta"
				SeqIO.write(sequences,fastafilename,"fasta")
				threshold,minthreshold=ComputeVariation(fastafilename,mincoverage)
				os.system("rm " + fastafilename)
			else:
				threshold,minthreshold=EvaluateVariation(taxonname,sequences,mincoverage)
			currentvariation={'Threshold': threshold, 'MinThreshold': minthreshold,'NumberOfSequences': len(sequences)}
			variations[taxonname]=currentvariation
		i=i+1	
	#write to file
	with open(variationfilename,"w") as json_file:
		json.dump(variations,json_file,encoding='latin1')
		
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

def ComputeBestBLASTscore(query,reference,mincoverage):
	indexed_query= IndexSequences(query)

	#load sequeces from the fasta files
	queryrecords = list(SeqIO.parse(indexed_query, "fasta"))
	#refrecords = list(SeqIO.parse(reference, "fasta"))

	bestscorelist =[0] * len(queryrecords)
	bestsimlist =[0] * len(queryrecords)
	bestcoveragelist =[0] * len(queryrecords)
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
			bestsimlist[i]=sim
			bestcoveragelist[i]=coverage
	os.system("rm " + indexed_query)		
	os.system("rm out.txt")
	return bestrefidlist,bestscorelist,bestsimlist,bestcoveragelist

def GetBestMatchLabels(trainclassification,trainseqIDs,bestmatchlist,bestscorelist,opthreshold,variation):
	bestlabels=[]
	i=0
	count=0
	for seqid in bestmatchlist:
		if  seqid in trainseqIDs:
			index=trainseqIDs.index(seqid)
			classname=unicode(trainclassification[index],'latin1')
			opt=0
			if optthreshold >= 0:
				opt=optthreshold
			else:
				opt=variation[classname]['Threshold']
#			if threshold < optthreshold:
#				opt=max(optthreshold -0.05,threshold)
			if bestscorelist[i] >= opt:
				bestlabels.append(trainclassification[index])
				count=count+1
			else:
			#no identification
				bestlabels.append("")
		else: #no identification
			bestlabels.append("")
		i=i+1
	return bestlabels,count

def SavePrediction(trainclassification,testclassification,testseqIDs,pred_labels,bestscorelist,bestsimlist,bestcoveragelist,bestrefidlist,opt,variation,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tClassification\tPrediction\tSimilarity core\tBLAST similarity score\tBLAST coverage\tBest match ID\tOptimal threshold\tThreshold\tMin threshold\tNumber of reference sequences\n")
	for i in range(0,len(testseqIDs)):
		testlabel=""
		if len(testclassification) > i:
			testlabel=testclassification[i]	
			predlabel=pred_labels[i]
			minthreshold=0
			threshold=0
			numberofsequences=0
			if predlabel!="":
				predlabel=unicode(predlabel,'latin1')
				if variation !={}:
					threshold=variation[predlabel]['Threshold']
					minthreshold=variation[predlabel]['MinThreshold']
					numberofsequences=variation[predlabel]['NumberOfSequences']
		output.write(str(i) + "\t" + str(testseqIDs[i]) + "\t" + testlabel + "\t"  + pred_labels[i] + "\t" + str(bestscorelist[i]) + "\t" + str(bestsimlist[i]) + "\t" + str(bestcoveragelist[i]) +  "\t" + bestrefidlist[i] + "\t" + str(opt) + "\t" + str(threshold) + "\t" + str(minthreshold) + "\t" + str(numberofsequences) + "\n")
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
trainclassification,classes,classnames,trainlevel= LoadClassification(trainseqIDs,trainseqrecords,classificationfilename, classificationposition)
testclassification,testclasses,testclassnames,testlevel= LoadClassification(testseqIDs,testseqrecords,classificationfilename, classificationposition)
		
#search for a best match of a test sequence in a train dataset
bestmatchlist,bestscorelist,bestsimlist,bestcoveragelist=ComputeBestBLASTscore(testdataset,traindataset,mincoverage)
#compute variation for reference sequences
variation={}
if optthreshold<0:
	if jsonvariationfilename=="":
		jsonvariationfilename = GetBase(traindataset) + "." + str(classificationposition) + ".variation"
		if not os.path.isfile(jsonvariationfilename):
			ComputeVariations(jsonvariationfilename,classes,classnames,mincoverage)
if jsonvariationfilename !=None: 
	if os.path.exists(jsonvariationfilename)==True:
		with open(jsonvariationfilename) as variation_file:
			variation = json.load(variation_file,encoding='latin1')	

#Get the best label for a test sequence based on its best match
bestlabels,count=GetBestMatchLabels(trainclassification,trainseqIDs,bestmatchlist,bestscorelist,optthreshold,variation)

#Save prediction by searching 
if reportfilename==None or reportfilename=="":	
	reportfilename=GetBase(testdataset) + "." + trainlevel + ".blast.classified"
SavePrediction(trainclassification,testclassification,testseqIDs,bestlabels,bestscorelist,bestsimlist,bestcoveragelist,bestmatchlist,optthreshold,variation,reportfilename)
print("Number of classified sequences: " + str(count))
print("The result is saved in file  " + reportfilename + ".")

