#!/usr/bin/env python
# FILE: classifyCNN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
from sklearn.model_selection import StratifiedKFold
#from sklearn.preprocessing import StandardScaler
import sys
import os
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Convolution1D
from keras.datasets import mnist
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution1D, MaxPooling1D
from keras.utils import np_utils
from keras import backend as K
from keras.models import load_model
import numpy as np
import json
from Bio import SeqIO
import random
import multiprocessing
nproc=multiprocessing.cpu_count()

modelname=sys.argv[1] #folder containing the classifier name and config file
testfastafilename=sys.argv[2]
minprobaforBlast=1.1 #search for best match of the sequences by BLast only when minprobabilityforBlast <=1.0
if len(sys.argv)>3:
	minprobaforBlast=float(sys.argv[3])
mincoverage=300 #for ITS sequences, used only for comparing the sequences with BLAST
if len(sys.argv)>4:
	mincoverage=int(sys.argv[4])
jsonvariationfilename="" #optional	
if len(sys.argv)>5:
	jsonvariationfilename=sys.argv[5]

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
	jsonvariationfilename=""
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
		if texts[0]=="Variation filename":
			jsonvariationfilename=texts[1].rstrip()	
	return jsonfilename,jsonvariationfilename,classifiername,classificationfilename,classificationpos,k,data_max

def loadData(matrixfilename,data_max,classificationfilename,classificationlevel):
	#load classification
	level=""
	seqids=[]
	records= open(classificationfilename)
	classification=[]
	for record in records:
		texts=record.split("\t")
		if record.startswith("#"):
			if classificationlevel < len(texts):
				level=texts[classificationlevel].rstrip()
			continue 
		seqid=texts[0].replace(">","").rstrip()
		classname=""
		if classificationlevel < len(texts):
			classname=texts[classificationlevel].rstrip()
		if classname !="":
			seqids.append(seqid)
			classification.append(classname)
	records.close()
	classificationset=set(classification)
	classes=list(classificationset)
	testvectors= list(open(matrixfilename, "r"))
	testvectors=testvectors[1:]
	X=[]
	testseqIDs=[]
	taxa=[]
	for vector in testvectors:
		elements=vector.split(",")
		seqid=elements[0]
		testseqIDs.append(seqid)
		taxonname=""
		if seqid in seqids:
			taxonname= classification[seqids.index(seqid)]
		taxa.append(taxonname)
		X.append(elements[1:])

	X=np.array(X,dtype=float)
	Y=[]
	i=0
	for taxonname in taxa:
		if taxonname in classes:
			Y.append(classes.index(taxonname))
		else:
			Y.append(-1)
		i=i+1
	Y=np.array(Y,dtype=int)	
	#data_max= np.amax(X)
	if data_max >0:
		X = X/data_max
	return testseqIDs,X,Y, classes,len(classes),len(X[0]),level	

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

def ComputeVariation(reffilename,mincoverage):
	#load sequeces from the fasta files
	records = list(SeqIO.parse(reffilename, "fasta"))
	scorematrix=ComputeBLASTscoreMatrix(reffilename,records,mincoverage)
	scorelist=[]
	for i in range(0,len(scorematrix)-2):
		for j in range(i+1,len(scorematrix)-1):
			if i!=j :
				scorelist.append(scorematrix[i][j])
	threshold=1
	minthreshold=1		
	if len(scorelist) >0:
		x=np.array(scorelist)
		minthreshold=np.min(x)
		threshold=np.median(x)
	return threshold,minthreshold

def EvaluateVariation(taxonname,sequences,mincoverage):
	thresholds=[]
	minthresholds=[] 
	for i in range(0,10):
		n=int(len(sequences)/10)
		selectedindexes=random.sample(range(0, len(sequences)), k=n)
		fastafilename=taxonname.replace(" ","_") + ".fasta"
		CreateFastaFileForSequences(fastafilename,sequences,selectedindexes)
		threshold,minthreshold=ComputeVariation(fastafilename,mincoverage)
		os.system("rm " + fastafilename)
		thresholds.append(threshold)
		minthresholds.append(minthreshold)
	threshold = np.median(np.array(thresholds))
	minthreshold =  np.median(np.array(minthresholds))
	return threshold,minthreshold

def ComputeVariations(variationfilename,classeswithsequences,mincoverage):
	#create json dict
	variations={}
	for taxonname in classeswithsequences.keys():
		sequences=classeswithsequences[taxonname]
		if len(sequences) >0:
			threshold=0
			minthreshold=0
			if len(sequences) < 100:				
				fastafilename=taxonname.replace(" ","_") + ".fasta"
				CreateFastaFileForSequences(fastafilename,sequences,[])
				threshold,minthreshold=ComputeVariation(fastafilename,mincoverage)
				os.system("rm " + fastafilename)
			else:
				threshold,minthreshold=EvaluateVariation(taxonname,sequences,mincoverage)
			currentvariation={'Threshold': threshold, 'MinThreshold': minthreshold,'NumberOfSequences': len(sequences)}
			variations[taxonname]=currentvariation
	#write to file
	with open(variationfilename,"w") as json_file:
		json.dump(variations,json_file,encoding='latin1')
		
def CreateFastaFileForSequences(fastafilename,sequences,selectedindexes):
	fastafile=open(fastafilename,"w")
	if len(selectedindexes)==0:
		for sequence in sequences:
			seq=sequence['seq']
			seqid=sequence['id']
			fastafile.write(">" + seqid + "\n")
			fastafile.write(seq + "\n")
	else:
		for index in selectedindexes:
			sequence=sequences[index]
			seq=sequence['seq']
			seqid=sequence['id']
			fastafile.write(">" + seqid + "\n")
			fastafile.write(seq + "\n")
	fastafile.close()

def CreateFastaFile(taxonname,classeswithsequences,seqno):
	taxonname=unicode(taxonname,errors='ignore')
	fastafilename=""
	sequences=[]
	if taxonname in classeswithsequences.keys():
		sequences=classeswithsequences[taxonname]
		if len(sequences) >0:
			fastafilename=taxonname.replace(" ","_").replace(".","_") + ".fasta"
			fastafile=open(fastafilename,"w")
			if (seqno >0) and (len(sequences) > seqno):
				#select randomly 100 sequences to compare
				selectedlist=random.sample(range(0, len(sequences)), k=seqno)
				for i in selectedlist:
					sequence=sequences[i]
					seqid=sequence['id']
					seq=sequence['seq']
					fastafile.write(">" + seqid + "\n")
					fastafile.write(seq + "\n")
			else:	
				for sequence in sequences:
					seqid=sequence['id']
					seq=sequence['seq']
					fastafile.write(">" + seqid + "\n")
					fastafile.write(seq + "\n")
			fastafile.close()
	return fastafilename,len(sequences)

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

def SavePrediction(minprobaforBlast,mincoverage,classes,classeswithsequences,variation,testseqIDs,testseqrecords,testlabels,pred_labels,probas,outputname):
	output=open(outputname,"w")
	if minprobaforBlast <=1.0:
		output.write("Index\tSequenceID\tClassification\tPrediction\tProbability\tlocal BLAST score\tlocal BLAST sim\tlocal BLAST coverage\tNumber of sequences to be compared\tBest match ID\tThreshold\tMin. threshold\n")
	else:
		output.write("Index\tSequenceID\tClassification\tPrediction\tProbability\n")
	i=0
	for seqid in testseqIDs:
		proba =probas[i][pred_labels[i]]
		giventaxonname=""
		if testlabels[i] >=0:
			giventaxonname = classes[testlabels[i]]
			giventaxonname=giventaxonname.encode('ascii', 'ignore')
		predictedname = classes[pred_labels[i]]
		if minprobaforBlast <= 1.0 and proba >= minprobaforBlast:
			reffilename,numberofsequences=CreateFastaFile(predictedname,classeswithsequences,0)
			bestrefid,bestlocalscore,bestlocalsim,bestlocalcoverage=ComputeBestLocalBLASTScore(testseqrecords[i],reffilename,mincoverage)
			predictedname=unicode(predictedname,'latin1')
			threshold=0
			minthreshold=0
			if predictedname in variation.keys():
					threshold=variation[predictedname]['Threshold']
					minthreshold=variation[predictedname]['MinThreshold']
			predictedname=predictedname.encode('ascii', 'ignore')
			output.write(str(i) + "\t" + seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t" + str(proba) + "\t" + str(bestlocalscore) + "\t" + str(bestlocalsim) + "\t" + str(bestlocalcoverage) + "\t" + str(numberofsequences) + "\t" + bestrefid + "\t" + str(threshold) + "\t" + str(minthreshold) + "\n")			
			os.system("rm " + reffilename)	
		else:
			output.write(str(i) + "\t" + seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t" + str(proba) + "\n")
		i=i+1
	output.close()
	
if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]
	#load config of the model
	jsonfilename,config_jsonvariationfilename,classifiername,classificationfilename,classificationlevel,k,data_max=LoadConfig(modelname)
	if jsonvariationfilename=="":
		jsonvariationfilename=config_jsonvariationfilename
	#represent sequences of the test dataset as k-mer vector
	testfilename=GetBase(testfastafilename)
	matrixfilename=testfilename + "." + str(k) + ".matrix"
	command=path + "fasta2matrix.py " +  str(k) + " " + testfastafilename + " " + matrixfilename
	os.system(command)
	testseqIDs,testinputs,testlabels,classes,nb_classes,input_length,level=loadData(matrixfilename,data_max,classificationfilename,classificationlevel)
	#load model
	model = load_model(classifiername)
	#predict labels for test dataset
	testinputs = testinputs.reshape(testinputs.shape + (1,))
	pred_labels = model.predict_classes(testinputs,verbose = 0)
	probas = model.predict_proba(testinputs)
	classeswithsequences={}
	variation={}
	testseqrecords=[]
	if minprobaforBlast <=1.0:
		testseqrecords=list(SeqIO.parse(testfastafilename, "fasta"))
		#load classes
		with open(jsonfilename) as json_file:
			   classeswithsequences = json.load(json_file)
#		#load variation
#		if not os.path.exists(jsonvariationfilename):
#			basename=modelname
#			if "/" in modelname:
#				basename=modelname[modelname.rindex("/")+1:]
#			jsonvariationfilename=modelname + "/" + basename + ".variation"
#			ComputeVariations(jsonvariationfilename,classeswithsequences,mincoverage)
#			#save into the contig file
#			configfilename=modelname + "/" + basename + ".config"
#			configfile=open(configfilename,"a")
#			configfile.write("Variation filename: " + jsonvariationfilename + "\n")
#			configfile.close()
		if os.path.exists(jsonvariationfilename):	   
			with open(jsonvariationfilename) as variation_file:
				variation = json.load(variation_file)
			   
	#save prediction
	basename=GetBase(classifiername)
	if "/" in basename:
		basename=basename[basename.rindex("/")+1:]
	reportfilename=GetBase(testfastafilename) +  "." + basename + ".classified"
	
	SavePrediction(minprobaforBlast,mincoverage,classes,classeswithsequences,variation,testseqIDs,testseqrecords,testlabels,pred_labels,probas,reportfilename)
	print("The result is saved in the file: " + reportfilename)
