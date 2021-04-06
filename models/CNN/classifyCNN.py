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
from keras.models import load_model
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
	classificationpos=0
	data_max=0
	k=6
	configfile=open(configfilename)
	for line in configfile:
		texts=line.split(": ")
		if texts[0]=="Classifier name":
			classifiername=texts[1].rstrip()
		if texts[0]=="K-mer number":
			k=int(texts[1].rstrip())
		if texts[0]=="Classification filename":
			classificationfilename=texts[1].rstrip()
		if texts[0]=="Data max":
			data_max=float(texts[1].rstrip())
		if texts[0]=="Column number to be classified":
			classificationpos=int(texts[1].rstrip())
		if texts[0]=="Classes filename":
			jsonfilename=texts[1].rstrip()
	return jsonfilename,classifiername,classificationfilename,classificationpos,k,data_max

def loadClassification(classificationfilename,classificationpos):
	#load classification
	allseqids=[]
	classificationfile= open(classificationfilename)
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
	output.write("SequenceID\tGiven classification\tPrediction\tFull classification of prediction\tProbability\n")
	i=0
	keys=classdict.keys()
	for seqid in testseqIDs:
		proba =probas[i][pred_labels[i]]
		giventaxonname=testtaxa[i]
		predictedname =keys[pred_labels[i]]
		#predictedname=predictedname.encode('ascii', 'ignore')	
		classification=classdict[predictedname]['classification']
		output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + str(proba) + "\n")
		i=i+1
	output.close()
	
if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]
	#load config of the model
	classesfilename,classifiername,classificationfilename,classificationpos,k,data_max=LoadConfig(modelname)
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
	model = load_model(classifiername)
	#predict labels for test dataset
	testinputs = testinputs.reshape(testinputs.shape + (1,))
	pred_labels = model.predict_classes(testinputs,verbose = 0)
	probas = model.predict_proba(testinputs)
	#save prediction
	basename=GetBase(classifiername)
	if "/" in basename:
		basename=basename[basename.rindex("/")+1:]
	if reportfilename==None or reportfilename=="":	
		reportfilename=GetBase(testfastafilename) +  "." + basename + ".classified"
	
	SavePrediction(classdict,testseqIDs,testtaxa,pred_labels,probas,reportfilename)
	print("The result is saved in the file: " + reportfilename)
