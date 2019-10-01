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


modelname=sys.argv[1] #folder containing the classifier name and config file
testfastafilename=sys.argv[2]
minprobaforBlast=-1 #search for best match of the sequences by BLast only when minprobabilityforBlast >=0
if len(sys.argv)>3:
	minprobaforBlast=float(sys.argv[3])

mincoverage=300 #for ITS sequences, used only for comparing the sequences with BLAST

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
		seqids.append(seqid)
		classification.append(texts[classificationlevel].rstrip())
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

def CreateFastaFile(taxonname,classeswithsequences):
	taxonname=unicode(taxonname,errors='ignore')
	fastafilename=""
	sequences=[]
	if taxonname in classeswithsequences.keys():
		sequences=classeswithsequences[taxonname]
		if len(sequences) >0:
			fastafilename=taxonname.replace(" ","") + ".fasta"
			fastafile=open(fastafilename,"w")
			if len(sequences) >100:
				#select randomly 500 sequences to compare
				selectedlist=random.sample(range(0, len(sequences)), k=100)
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

def GetBestLocalSim(localsims,localcoverages,mincoverage):
	i=0
	bestlocalsim=0
	bestlocalcoverage=0
	for sim in localsims:
		if localcoverages[i] >=mincoverage:
			if sim>bestlocalsim:
				bestlocalsim=sim
				bestlocalcoverage=localcoverages[i]
		i=i+1		
	if bestlocalsim==0 and len(localcoverages) > 0:
		bestlocalcoverage=max(localcoverages)
		bestlocalsim=localsims[localcoverages.index(bestlocalcoverage)]
	return bestlocalsim,bestlocalcoverage

def ComputeBLASTScore(testrecord,predictedname,classeswithsequences,mincoverage):
	#Create fasta file of the test record 
	queryname=testrecord.id + "." + "fasta"
	if "|" in testrecord.id:
		queryname=testrecord.id.split("|")[0] + "." + "fasta"
	SeqIO.write(testrecord, queryname , "fasta")
	#Create fasta file of predictedname
	refname,numberofsequences=CreateFastaFile(predictedname,classeswithsequences)
	if refname=="":
		return "",0,0,0,0,0,0
	#Blast the test record to fastafilename:
	#makedbcommand = "makeblastdb -in " + refname + " -dbtype \'nucl\' " +  " -out db"
	#os.system(makedbcommand)
	blastoutput = testrecord.id 
	if "|" in testrecord.id:
		blastoutput=blastoutput.split("|")[0]
	blastoutput= blastoutput + ".blast.out"
	#blastcommand = "blastn -query " + queryname + " -db  db -outfmt 6 -out " + blastoutput
	blastcommand = "blastn -query " + queryname + " -subject " + refname  + " -outfmt 6 -out " + blastoutput
	os.system(blastcommand)
	#read output of Blast
	blastoutputfile = open(blastoutput)
	bestscore=0
	bestsim=0
	bestcoverage =0
	bestlocalsim=0
	bestlocalcoverage=0
	bestrefid=""
	refid = ""
	queryid=""
	minposition = 0
	maxposition = 0
	refminposition = 0
	refmaxposition = 0
	positions  = []
	identities =[]
	localsims=[]
	localcoverages=[]
	score=0
	identity =0
	for line in blastoutputfile:
		words = line.split("\t")
		currentqueryid = words[0]
		currentrefid = words[1]
		if (refid != currentrefid) or (currentqueryid !=queryid):
			if (refid !="") and (queryid != ""):
				coverage = max(abs(maxposition - minposition) + 1,abs(refmaxposition -refminposition) + 1) + min(minposition,refminposition)-1 
				#coverage = max(abs(maxposition - minposition) + 1,abs(refmaxposition -refminposition) + 1) 
				#compute identity
				identity =0 
				for pos in range(minposition,maxposition + 1):
					if pos in positions:
						indexofpos = positions.index(pos)
						identitiesatpos=identities[indexofpos]
						identity = identity + max(identitiesatpos)/100
				sim=0
				if coverage > 0:
					sim  = float(identity) / float(coverage)
				score = sim
				localsim,localcoverage=GetBestLocalSim(localsims,localcoverages,mincoverage)
				if localcoverage>=mincoverage:
					score=localsim	
				elif coverage < mincoverage:
					score = float(score * coverage) /mincoverage
				if bestscore < score:
					bestscore= score
					bestsim= sim
					bestcoverage= coverage
					bestrefid=refid
					bestlocalsim=localsim
					bestlocalcoverage=localcoverage
			queryid=currentqueryid
			refid = currentrefid
			minposition = 0
			maxposition =0
			refminposition = 0
			refmaxposition =0 
			refid = currentrefid
			identity = 0
			positions = []
			identities =[]
			localsims=[]
			localcoverages=[]
		refpos1 =  int(words[8])
		refpos2 = int(words[9])
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		localsims.append(float(iden)/100)
		localcoverages.append(abs(pos2-pos1))
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
	coverage = max(abs(maxposition - minposition) + 1,abs(refmaxposition -refminposition) + 1) + min(minposition,refminposition)-1 
	#coverage = max(abs(maxposition - minposition) + 1,abs(refmaxposition -refminposition) + 1)
	#compute identity
	identity =0 
	for pos in range(minposition,maxposition + 1):
		if pos in positions:
			indexofpos = positions.index(pos)
			identitiesatpos=identities[indexofpos]
			identity = identity + max(identitiesatpos)/100
	if coverage > 0:
		sim  = identity / coverage
		score = sim
		localsim,localcoverage=GetBestLocalSim(localsims,localcoverages,mincoverage)
		if localcoverage>=mincoverage:
			score=localsim	
		elif coverage < mincoverage:
			score = float(score * coverage) /mincoverage
		if bestscore < score:
			bestscore= score
			bestsim= sim
			bestcoverage= coverage
			bestrefid=refid
			bestlocalsim=localsim
			bestlocalcoverage=localcoverage
	os.system("rm " + blastoutput)
	os.system("rm " + queryname)
	os.system("rm " + refname)
	return bestrefid,bestscore,bestsim,bestcoverage,bestlocalsim,bestlocalcoverage,numberofsequences

def SavePrediction(minprobaforBlast,mincoverage,classes,classeswithsequences,testseqIDs,testseqrecords,testlabels,pred_labels,probas,outputname):
	output=open(outputname,"w")
	if minprobaforBlast >= 0:
		output.write("Index\tSequenceID\tClassification\tPrediction\tProbability\tBLAST score\tBLAST sim\tBLAST coverage\tBLAST local sim\tBLAST local coverage\tNumber of sequences to be compared\tBest match ID\n")
	else:
		output.write("Index\tSequenceID\tClassification\tPrediction\tProbability\t\n")	
	acc=0
	i=0
	for seqid in testseqIDs:
		proba =probas[i][pred_labels[i]]
		giventaxonname=""
		if testlabels[i] >=0:
			giventaxonname = classes[testlabels[i]]
		predictedname = classes[pred_labels[i]]
		if minprobaforBlast >= 0:
			bestrefseqID,score,sim,coverage,localsim,localcoverage,numberofsequences=ComputeBLASTScore(testseqrecords[i],predictedname,classeswithsequences,mincoverage)
			output.write(str(i) + "\t" + seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t" + str(proba) + "\t" + str(score) + "\t" + str(sim) + "\t" + str(coverage) + "\t" + str(localsim) + "\t" + str(localcoverage) + "\t" + str(numberofsequences) + "\t" + bestrefseqID + "\n")			
		else:
			output.write(str(i) + "\t" + seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t" + str(proba) + "\n")
		i=i+1
	output.close()
	
if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]
	#load config of the model
	jsonfilename,classifiername,classificationfilename,classificationlevel,k,data_max=LoadConfig(modelname)
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
	testseqrecords=[]
	if minprobaforBlast >=0:
		#load classes
		with open(jsonfilename) as json_file:
			   classeswithsequences = json.load(json_file)
		testseqrecords=list(SeqIO.parse(testfastafilename, "fasta"))
	#save prediction
	basename=GetBase(classifiername)
	if "/" in basename:
		basename=basename[basename.rindex("/")+1:]
	reportfilename=GetBase(testfastafilename) +  "." + basename + ".classified"
	SavePrediction(minprobaforBlast,mincoverage,classes,classeswithsequences,testseqIDs,testseqrecords,testlabels,pred_labels,probas,reportfilename)
	print("The result is saved in the file: " + reportfilename)
