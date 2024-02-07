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
import tensorflow as tf
import numpy as np
import json
#from Bio import SeqIO
#import random
import multiprocessing
parser=argparse.ArgumentParser(prog='classifyCNN.py',  
							   usage="%(prog)s [options] -i fastafile -c classifiername -mp minproba -mc mincoverage -j variationjsonfilename",
							   description='''Script that classifies the sequences of the fasta files using CNN model. The classified sequences with a probability less than minproba will be verified using BLAST. The mincoverage is given for BLAST comparision. The json file is the the variation of the sequences within each group of the training dataset. This file is created during the training of the model, and used optionally for the verification of the classification. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The folder name containing the model and associated files.') #optional
parser.add_argument('-c','--classifier', required=True, help='the folder containing the classifier.')

args=parser.parse_args()
testfastafilename= args.input
modelname=args.classifier
reportfilename=args.out

nproc=multiprocessing.cpu_count()

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadConfig(modelname):
	if modelname[len(modelname)-1]=="/":
		modelname=modelname[:-1]
	basename=modelname
	if "/" in modelname:
		basename=modelname[modelname.rindex("/")+1:]
	configfilename=modelname + "/" + basename + ".config"
	classifiername=""
	jsonfilename=""
	classificationfilename=""
	data_max=0
	k=6
	rank=""
	classificationpos=-1
	configfile=open(configfilename)
	for line in configfile:
		texts=line.split(": ")
		if texts[0]=="Classifier name":
			classifiername=texts[1].rstrip()
			if "/" in classifiername:
				classifiername=modelname + "/" + classifiername[classifiername.rindex("/")+1:]
		if texts[0]=="K-mer number":
			k=int(texts[1].rstrip())
		if texts[0]=="Classification filename":
			classificationfilename= texts[1].rstrip()
			if "/" in classificationfilename:
				classificationfilename=modelname + "/" + classificationfilename[classificationfilename.rindex("/")+1:]
		if texts[0]=="Data max":
			data_max=float(texts[1].rstrip())
		if texts[0]=="Rank to be classified":
			rank=texts[1].rstrip()
		if texts[0]=="Column number to be classified":
			classificationpos=int(texts[1].rstrip())	
		if texts[0]=="Classes filename":
			jsonfilename= texts[1].rstrip()
			if "/" in jsonfilename:
				jsonfilename=modelname + "/" + jsonfilename[jsonfilename.rindex("/")+1:]
	if classificationpos==-1 and rank!="":
		classificationpos=GetClassificationPos(classificationfilename,rank)		
	return jsonfilename,classifiername,classificationfilename,rank,classificationpos,k,data_max

def loadClassification(classificationfilename,classificationpos):
	#load classification
	allseqids=[]
	classificationfile= open(classificationfilename, errors='ignore')
	classifications=[]
	for line in classificationfile:
		texts=line.split("\t")
		if line.startswith("#"):
			continue 
		seqid=texts[0].replace(">","").rstrip()
		classname=""
		if classificationpos < len(texts):
			classname=texts[classificationpos].rstrip()
		if classname !="":
			allseqids.append(seqid)
			classifications.append(classname)
	classificationfile.close()
	return classifications,allseqids

def loadData(matrixfilename,data_max,classifications,allseqids):
	#load vectors from the matrix
	testvectors= list(open(matrixfilename, "r"))
	testvectors=testvectors[1:]
	X=[]
	testseqIDs=[]
	testtaxa=[]
	for vector in testvectors:
		elements=vector.split(",")
		seqid=elements[0]
		testseqIDs.append(seqid)
		taxonname=""
		if seqid in allseqids:
			taxonname= classifications[allseqids.index(seqid)]
		testtaxa.append(taxonname)
		X.append(elements[1:])
	X=np.array(X,dtype=float)
	if data_max >0:
		X = X/data_max
	return testseqIDs,X,testtaxa

def SavePrediction(classdict,testseqIDs,testtaxa,pred_labels,probas,outputname):
	output=open(outputname,"w")
	output.write("#SequenceID\tGiven label\tPrediction\tFull classification\tProbability\n")
	i=0
	keys=list(classdict.keys())
	for seqid in testseqIDs:
		proba =probas[i][pred_labels[i]]
		giventaxonname=testtaxa[i]
		predictedname =keys[pred_labels[i]]
		classification=classdict[predictedname]['classification']
		output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + str(proba) + "\n")
		i=i+1
	output.close()

def GetClassificationPos(classificationfilename,rank):
	classificationfile=open(classificationfilename, errors='ignore')
	header=next(classificationfile)
	texts=header.split("\t")
	i=0
	classificationpos=-1
	for text in texts:
		text=text.rstrip()
		if text.lower()==rank.lower():
			classificationpos=i
		i=i+1	
	classificationfile.close()
	return classificationpos

if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]
	#load config of the model
	classesfilename,classifiername,classificationfilename,rank,classificationpos,k,data_max=LoadConfig(modelname)
	#load ref class dict
	classdict={}
	with open(classesfilename) as classesfile:
		classdict = json.load(classesfile)
	#load classification
	classifications,allseqids=loadClassification(classificationfilename,classificationpos)	
	#represent sequences of the test dataset as k-mer vector
	testfilename=GetBase(testfastafilename)
	matrixfilename=testfilename + "." + str(k) + ".matrix"
	command=path + "fasta2matrix.py " +  str(k) + " " + testfastafilename + " " + matrixfilename
	os.system(command)
	testseqIDs,testinputs,testtaxa=loadData(matrixfilename,data_max,classifications,allseqids)
	#load model
	model = tf.keras.models.load_model(classifiername)
	# Predict labels for test dataset
	testinputs = testinputs.reshape(testinputs.shape + (1,))
	probas = model.predict(testinputs, verbose=0)
	pred_labels = np.argmax(probas, axis=1)
	#save prediction
	basename=GetBase(classifiername)
	if "/" in basename:
		basename=basename[basename.rindex("/")+1:]
	if reportfilename==None or reportfilename=="":	
		reportfilename=GetBase(testfastafilename) +  "." + basename + ".out"
		
	SavePrediction(classdict,testseqIDs,testtaxa,pred_labels,probas,reportfilename)
	print("The result is saved in the file: " + reportfilename)
