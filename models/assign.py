#!/usr/bin/env python
# FILE: classifyCNN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
#from sklearn.model_selection import StratifiedKFold
#from sklearn.preprocessing import StandardScaler
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
#from keras.models import Sequential
#from keras.layers import Dense
#from keras.layers import Convolution1D
#from keras.datasets import mnist
#from keras.layers import Dropout, Activation, Flatten
#from keras.layers import Convolution1D, MaxPooling1D
#from keras.utils import np_utils
#from keras import backend as K
#from keras.models import load_model
#import numpy as np
import json
from Bio import SeqIO
import random
import multiprocessing
parser=argparse.ArgumentParser(prog='assign.py',  
							   usage="%(prog)s [options] -i fastafile -p prediction -r referencefastafile -c classificationfile -mp minproba -mc mincoverage -cutoffs cutoffsfile -o output",
							   description='''Script that assigns the classified sequences of the prediction file to their BLAST best match based on the given cutoffs.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-r','--reference', required=True, help='the reference fasta file')
parser.add_argument('-o','--out', help='The assignment of the sequences.')
parser.add_argument('-p','--predictionfile', required=True, help='the prediction/classification of the sequences.')
parser.add_argument('-c','--classificationfile', required=True, help='the classification file in tab. format.')
parser.add_argument('-mp','--minproba', type=float, default=0, help='The minimum probability for verifying the classification results.')
parser.add_argument('-mc','--mincoverage', type=int, default=300, help='Minimum coverage required for the identitiy of the BLAST comparison.')
parser.add_argument('-m','--maxseqno', type=int, default=0, help='Maximum number of the sequences of the predicted taxon name from the classification file will be selected for the comparison to find the best match. If it is not given, all the sequences will be selected.')
parser.add_argument('-cutoffs','--cutoffs', help='The json file containing the cutoffs of the predicted taxa.')

args=parser.parse_args()
fastafilename= args.input
referencefastafilename= args.reference
outputname=args.out
minprobaforBlast=args.minproba
mincoverage = args.mincoverage
cutoffsfilename=args.cutoffs
classificationfilename=args.classificationfile
predictionfilename=args.predictionfile
maxseqno=args.maxseqno


nproc=multiprocessing.cpu_count()

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

def loadData(reffastafilename,classificationfilename,classificationpos):
	#load referencesequences
	refdict=SeqIO.to_dict(SeqIO.parse(reffastafilename, "fasta"))
	#load classification
	classificationfile= open(classificationfilename)
	classeswithsequences={}
	classifications={}
	header=""
	for line in classificationfile:
		texts=line.rstrip().split("\t")
		if line.startswith("#"): 
			header=line.rstrip()					 
			continue				 
		if classificationpos < len(texts):
			classname=texts[classificationpos]
			classification=GetTaxonomicClassification(classificationpos,header,texts)
			seq=refdict[texts[0]]
			if classname in classeswithsequences.keys():
				classeswithsequences[classname].append(seq)
			else:
				classeswithsequences.setdefault(classname,[seq])
			if classname in classifications.keys():
				if len(classifications[classname]) < len(classification):
					classifications[classname]=classification
			else:
				classifications.setdefault(classname,classification)	
	return classeswithsequences,classifications

def CreateFastaFile(taxonname,classeswithsequences,maxseqno):
	if sys.version_info[0] < 3:
		taxonname=unicode(taxonname,errors='ignore')
	newfastafilename=""
	sequences=[]
	if taxonname in classeswithsequences.keys():
		sequences=classeswithsequences[taxonname]
		if len(sequences) >0:
			newfastafilename=taxonname.replace(" ","_").replace(".","_") + ".fasta"
			fastafile=open(newfastafilename,"w")
			if (maxseqno >0) and (len(sequences) > maxseqno):
				#select randomly 100 sequences to compare
				selectedlist=random.sample(range(0, len(sequences)), k=maxseqno)
				for i in selectedlist:
					sequence=sequences[i]
					seqid=sequence.id
					seq=str(sequence.seq)
					fastafile.write(">" + seqid + "\n")
					fastafile.write(seq + "\n")
			else:	
				for sequence in sequences:
					seqid=sequence.id
					seq=str(sequence.seq)
					fastafile.write(">" + seqid + "\n")
					fastafile.write(seq + "\n")
			fastafile.close()
	return newfastafilename,len(sequences)

def ComputeBestLocalBLASTScore(testrecord,reffilename,mincoverage):
	#Create fasta file of the test record 
	queryname=testrecord.id + "." + "fasta"
	if "|" in testrecord.id:
		queryname=testrecord.id.split("|")[0] + "." + "fasta"
	SeqIO.write(testrecord, queryname , "fasta")
	#Create fasta file of predictedname
	if reffilename=="":
		return "",0,0,0
	#Blast the test record to fastafilename:
	#makedbcommand = "makeblastdb -in " + refname + " -dbtype \'nucl\' " +  " -out db"
	#os.system(makedbcommand)
	blastoutput = testrecord.id 
	if "|" in testrecord.id:
		blastoutput=blastoutput.split("|")[0]
	blastoutput= blastoutput + ".blast.out"
	#blastcommand = "blastn -query " + queryname + " -db  db -outfmt 6 -out " + blastoutput
	blastcommand = "blastn -query " + queryname + " -subject " + reffilename  + " -task blastn-short -outfmt 6 -out " + blastoutput
	os.system(blastcommand)
	#read output of Blast
	blastoutputfile = open(blastoutput)
	bestlocalscore=0
	bestlocalsim=0
	bestlocalcoverage=0
	bestrefid=""
	for line in blastoutputfile:
		words = line.split("\t")
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		refid=words[1]
		score=sim
		if coverage < mincoverage:
				score=float(score * coverage)/mincoverage
		if score > bestlocalscore:
			bestrefid=refid
			bestlocalscore=score
			bestlocalsim=sim
			bestlocalcoverage=coverage
	os.system("rm " + blastoutput)
	os.system("rm " + queryname)
	return bestrefid,bestlocalscore,bestlocalsim,bestlocalcoverage

def Assign(classeswithsequences,classifications,classificationpos,rank,minprobaforBlast,mincoverage,cutoffs,seqids,seqdict,labels,pred_labels,probas,maxseqno,outputname):
	output=open(outputname,"w")
	output.write("Index\tSequenceID\tGiven classification\tPrediction\tFull classification of prediction\tProbability\tAssigned\tReferenceID\tCut-off\tConfidence\tBLAST score\tBLAST sim\tBLAST coverage\n")
	i=0
	for seqid in seqids:
		proba =probas[i]
		giventaxonname=labels[i]
		predictedname=pred_labels[i]
		assignment=""
		score=0
		sim=0
		coverage=0
		if proba >= minprobaforBlast:
			reffilename,numberofsequences=CreateFastaFile(predictedname,classeswithsequences,maxseqno)
			refid,score,sim,coverage=ComputeBestLocalBLASTScore(seqdict[seqid],reffilename,mincoverage)
			os.system("rm " + reffilename)	
			if sys.version_info[0] < 3:
				predictedname=unicode(predictedname,'latin1')
			cutoff=0
			confidence=0
			rank_cutoff=0
			rank_confidence=0
			if rank in cutoffs.keys():
				rank_cutoff=cutoffs[rank]['cut-off']
				rank_confidence=cutoffs[rank]['confidence']
			elif str(classificationpos) in cutoffs.keys():
				rank_cutoff=cutoffs[rank]['cut-off']
				rank_confidence=cutoffs[rank]['confidence']
			if predictedname in cutoffs.keys():
				cutoff=cutoffs[predictedname]['cut-off']
				confidence=cutoffs[predictedname]['confidence']
			if confidence < rank_confidence:
				cutoff=rank_cutoff
				confidence=rank_confidence
			if sys.version_info[0] < 3:	
				predictedname=predictedname.encode('ascii', 'ignore')
				classification=classifications[predictedname]
			assignment="no"		
			if score >=cutoff:
				assignment="yes"	
		output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + str(proba) + "\t" + assignment + "\t" + refid + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + str(score) + "\t" + str(sim) + "\t" + str(coverage) + "\n")			
		i=i+1
	output.close()
	
def LoadPrediction(predictionfilename):
	seqids=[]
	labels=[]
	pre_labels=[]
	classifications=[]
	probas=[]
	predictionfile=open(predictionfilename)
	next(predictionfile)
	for line in predictionfile:
		texts=line.rstrip().split("\t")
		seqids.append(texts[0])
		labels.append(texts[1])
		pre_labels.append(texts[2])
		classifications.append(texts[3])
		probas.append(float(texts[4]))
	return seqids,labels,pre_labels,classifications,probas
	
def GetClassificationpos(pred_labels,classificationfilename):
	classificationpos=0
	rank=""
	classificationfile=open(classificationfilename)
	for line in classificationfile:
		if line.startswith("#"):
			continue
		texts=line.rstrip().split("\t")	 
		i=0
		for text in texts:
			if text!="":
				if text in pred_labels:
					classificationpos=i
					break
			i=i+1	
		if classificationpos >0:
			break
	classificationfile.close()	
	
	return classificationpos,rank
if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]

	basename=GetBase(referencefastafilename)
	if "/" in basename:
		basename=basename[basename.rindex("/")+1:]
	if outputname==None or outputname=="":	
		outputname=GetBase(fastafilename) +  "." + basename + ".assigned"
	#load prediction
	seqids,labels,pred_labels,classifications,probas=LoadPrediction(predictionfilename)	
	#load reference sequences
	classificationpos,rank=GetClassificationpos(pred_labels,classificationfilename)
	classeswithsequences,classifications=loadData(referencefastafilename,classificationfilename,classificationpos)	
	#load seqdict	
	seqdict=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	cutoffs={}
	if cutoffsfilename!="" and cutoffsfilename!=None:
		with open(cutoffsfilename) as cutoffsfile:
			cutoffs = json.load(cutoffsfile)
	Assign(classeswithsequences,classifications,classificationpos,rank,minprobaforBlast,mincoverage,cutoffs,seqids,seqdict,labels,pred_labels,probas,maxseqno,outputname)
	print("The result is saved in the file: " + outputname)
