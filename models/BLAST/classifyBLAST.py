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
parser.add_argument('-t','--threshold', required=True, type=float, default=0.97, help='The threshold for the classification.')
parser.add_argument('-mc','--mincoverage', type=int, default=300, help='Optinal. Minimum coverage required for the identitiy of the BLAST comparison.')
parser.add_argument('-o','--out', help='The folder name containing the model and associated files.') #optional
parser.add_argument('-c','--classification',help='the classification file in tab. format.')#optinal
parser.add_argument('-p','--classificationpos', type=int, default=0, help='the classification position to load the classification.')#optional

args=parser.parse_args()
testdataset= args.input
traindataset = args.reference
optthreshold=args.threshold
mincoverage = args.mincoverage
classificationfilename=args.classification
classificationposition=args.classificationpos
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

def GetTaxonomicClassification(classificationpos,header,texts):
	classification=""
	p_s=len(texts)
	p_g=len(texts)
	p_f=len(texts)
	p_o=len(texts)
	p_c=len(texts)
	p_p=len(texts)
	p_k=len(texts)
	i=0
	for text in header.split("\t"):
		text=text.rstrip()
		if text.lower()=="species":
			p_s=i
		elif text.lower()=="genus":
			p_g=i	
		elif text.lower()=="family":
			p_f=i	
		elif text.lower()=="order":
			p_o=i	
		elif text.lower()=="class":
			p_c=i	
		elif text.lower()=="phylum":
			p_p=i	
		elif text.lower()=="kingdom":
			p_k=i	
		i=i+1 
	species="s__"
	genus="g__"
	family="f__"
	order="o__" 
	bioclass="c__"
	phylum="p__"
	kingdom="k__"
	if p_s< len(texts):
		species="s__" + texts[p_s].rstrip()
	if p_g< len(texts):
		genus="g__" + texts[p_g].rstrip()	
	if p_f< len(texts):
		family="f__" + texts[p_f].rstrip()	
	if p_o< len(texts):
		order="o__" + texts[p_o].rstrip()
	if p_c< len(texts):
		bioclass="c__" + texts[p_c].rstrip()	
	if p_p< len(texts):
		phylum="p__" + texts[p_p].rstrip()
	if p_k< len(texts):
		kingdom="k__" + texts[p_k].rstrip()	
	if classificationpos==p_s:
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family + ";"+ genus+";"+species
	elif classificationpos==p_g:
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family + ";"+ genus
	elif classificationpos==p_f:
		classification=kingdom +";"+phylum + ";"+bioclass +";"+order+";"+family 
	elif classificationpos==p_o:
		classification=kingdom +";"+phylum + ";"+bioclass + ";"+order
	elif classificationpos==p_c:
		classification=kingdom +";"+phylum + ";"+bioclass 
	elif classificationpos==p_p:
		classification=kingdom +";"+phylum
	elif classificationpos==p_k:
		classification=kingdom
	else:
		classification=texts[classificationpos]
	return classification

def LoadClassification(seqIDs,seqrecords,classificationfilename,pos):
	classification=[""]*len(seqIDs)
	fullclassification=[""]*len(seqIDs)
	classes=[]
	classnames=[]
	level=str(pos)
	if classificationfilename == "":
		return classification
	classificationfile= list(open(classificationfilename, "r"))
	header=""
	for line in classificationfile:
		texts=line.split("\t")
		if line.startswith("#"):
			header=line		 
			level=texts[pos]
			continue 
		seqid = texts[0].replace(">","").rstrip()
		if seqid in seqIDs:
			index=seqIDs.index(seqid)
			classname=texts[pos].rstrip()
			classification[index]=classname
			fullclassification[index]=GetTaxonomicClassification(pos,header,texts)
			if classname in classnames:
				classid=classnames.index(classname)
				classes[classid].append(seqrecords[index])
			else:
				classnames.append(classname)
				seqs=[]
				seqs.append(seqrecords[index])
				classes.append(seqs)
	return fullclassification,classification,classes,classnames,level

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

def GetBestMatchLabels(trainclassification,trainseqIDs,bestmatchlist,bestscorelist,opthreshold):
	bestlabels=[]
	i=0
	count=0
	for seqid in bestmatchlist:
		if  seqid in trainseqIDs:
			index=trainseqIDs.index(seqid)
			classname=trainclassification[index]
			if sys.version_info[0] < 3:
				classname=unicode(trainclassification[index],'latin1')
			if bestscorelist[i] >= opthreshold:
				bestlabels.append(classname)
				count=count+1
			else:
			#no identification
				bestlabels.append("")
		else: #no identification
			bestlabels.append("")
		i=i+1
	return bestlabels,count

def SavePrediction(fulltrainclassification,trainclassification,testclassification,testseqIDs,pred_labels,bestscorelist,bestsimlist,bestcoveragelist,bestrefidlist,opt,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tGiven classification\tPrediction\tFull classification of prediction\tBLAST score\tBLAST similarity\tBLAST coverage\tBest match ID\tOptimal threshold\n")
	for i in range(0,len(testseqIDs)):
		testlabel=""
		if len(testclassification) > i:
			testlabel=testclassification[i]	
			predlabel=pred_labels[i]
			if predlabel!="" and sys.version_info[0] < 3:
				predlabel=unicode(predlabel,'latin1')
			classification=fulltrainclassification[trainclassification.index(predlabel)]
		output.write(str(i) + "\t" + str(testseqIDs[i]) + "\t" + testlabel + "\t"  + pred_labels[i] + "\t" + classification +"\t" + str(bestscorelist[i]) + "\t" + str(bestsimlist[i]) + "\t" + str(bestcoveragelist[i]) +  "\t" + bestrefidlist[i] + "\t" + str(opt) +"\n")
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
fulltrainclassification,trainclassification,classes,classnames,trainlevel= LoadClassification(trainseqIDs,trainseqrecords,classificationfilename, classificationposition)
fulltestclassification,testclassification,testclasses,testclassnames,testlevel= LoadClassification(testseqIDs,testseqrecords,classificationfilename, classificationposition)
		
#search for a best match of a test sequence in a train dataset
bestmatchlist,bestscorelist,bestsimlist,bestcoveragelist=ComputeBestBLASTscore(testdataset,traindataset,mincoverage)

#Get the best label for a test sequence based on its best match
bestlabels,count=GetBestMatchLabels(trainclassification,trainseqIDs,bestmatchlist,bestscorelist,optthreshold)

#Save prediction by searching 
if reportfilename==None or reportfilename=="":	
	reportfilename=GetBase(testdataset) + "." + trainlevel + ".blast.classified"
SavePrediction(fulltrainclassification,trainclassification,testclassification,testseqIDs,bestlabels,bestscorelist,bestsimlist,bestcoveragelist,bestmatchlist,optthreshold,reportfilename)
print("Number of classified sequences: " + str(count))
print("The result is saved in file  " + reportfilename + ".")

