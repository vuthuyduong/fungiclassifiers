#!/usr/bin/env python
# FILE: trainCNN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
from sklearn.model_selection import StratifiedKFold
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
import numpy as np
import json

fastafilename=sys.argv[1]
classificationfilename=sys.argv[2] #taxonomy file 
classificationlevel=int(sys.argv[3]) #the level of classification to get taxa from the taxonomy file
k = 6
if len(sys.argv) >4:
	k= int(sys.argv[4])

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]
		
def load_data(fastafilename,classificationfilename,classificationlevel):
	#load classification
	allseqids=[]
	records= open(classificationfilename)
	classification=[]
	level=""
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
			allseqids.append(seqid)
			classification.append(classname)
	records.close()
	classificationset=set(classification)
	classes=list(classificationset)
	#load fastafile, save a new fasta file containing only sequences having a classification
	fastafile=open(fastafilename)
	newfastafilename=GetBase(fastafilename) + "." + str(classificationlevel) + ".fasta"
	newfastafile=open(newfastafilename,"w")
	writetofile=False
	seq=""
	sequences=[]
	seqids=[]
	taxa=[]
	for line in fastafile:
		if line.startswith(">"):
			if writetofile==True:
				sequences.append(seq)
			writetofile=False
			seq=""
			seqid=line.split("|")[0].replace(">","").rstrip()		
			if seqid in allseqids:
				index=allseqids.index(seqid)
				taxonname=classification[index]
				if taxonname !="":
					#write to file
					newfastafile.write(">" + seqid + "\n")
					seqids.append(seqid)
					taxa.append(taxonname)
					writetofile=True
		else:
			if writetofile==True:
				newfastafile.write(line)
				seq=seq+line.rstrip()
	if writetofile==True:
		sequences.append(seq)
	fastafile.close()
	newfastafile.close()
	return newfastafilename, seqids, sequences, taxa, classes,level

def load_matrix(matrixfilename,seqids,sequences,taxa,classes):
	#load matrix
	vectors= list(open(matrixfilename, "r"))
	vectors=vectors[1:]
	X=[]
	Y=[]
	seqIDList=[]
	seqList=[]
	for i in range(0,len(classes)):
		seqIDList.append([])
		seqList.append([])
	for vector in vectors:
		elements=vector.split(",")
		seqid=elements[0]
		taxonname= taxa[seqids.index(seqid)]
		seq=sequences[seqids.index(seqid)]
		index=classes.index(taxonname)
		Y.append(index)
		X.append(elements[1:])
		if seqid not in seqIDList[index]:
			seqIDList[index].append(seqid)
			seqList[index].append(seq)
	X=np.array(X,dtype=float)
	Y=np.array(Y,dtype=int)
	data_max=0
	#data_max= np.amax(X)
	#X = X/data_max
	return X,Y,len(classes),len(X[0]),data_max,seqIDList,seqList

def create_model(nb_classes,input_length):
	model = Sequential()
	model.add(Convolution1D(5,5, border_mode='valid', input_dim=1,input_length=input_length)) #input_dim
	model.add(Activation('relu'))
	model.add(MaxPooling1D(pool_length=2,border_mode='valid'))
	model.add(Convolution1D(10, 5,border_mode='valid'))
	model.add(Activation('relu'))
	model.add(MaxPooling1D(pool_length=2,border_mode='valid'))
	model.add(Flatten())
	##
	##MLP
	model.add(Dense(500))
	model.add(Activation('relu'))
	model.add(Dropout(0.5))
	model.add(Dense(nb_classes))
	model.add(Activation('softmax'))
	model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])
	return model

def SaveConfig(configfilename,classifiername,fastafilename,jsonfilename,classificationfilename,classificationpos,kmer,data_max):
	if not classifiername.startswith("/"):
		classifiername=os.getcwd() + "/" + classifiername
	if not fastafilename.startswith("/"):
		fastafilename=os.getcwd() + "/" + fastafilename
	if not jsonfilename.startswith("/"):
		jsonfilename=os.getcwd() + "/" + jsonfilename
	if not classificationfilename.startswith("/"):
		classificationfilename=os.getcwd() + "/" + classificationfilename
	model="cnn"
	#save the config:classifierfilename, model, classificationfilename,classificationpos,k-mer
	configfile=open(configfilename,"w")
	configfile.write("Classifier name: " + classifiername + "\n")
	configfile.write("Model: " + model + "\n")
	configfile.write("Data max: " + str(data_max) + "\n")
	configfile.write("K-mer number: " + str(kmer) + "\n")
	configfile.write("Fasta filename: " + fastafilename + "\n")
	configfile.write("Classification filename: " + classificationfilename + "\n")
	configfile.write("Column number to be classified: " + str(classificationpos) + "\n")
	configfile.write("Classes filename: " + jsonfilename + "\n")
	configfile.close()	

def SaveClasses(jsonfilename,classnames,seqIDList,seqList):
	#create json dict
	taxadict={}
	i=0
	count=0
	for classname in classnames:	
		classname=unicode(classname,errors='ignore')
		currentclass=[]
		j=0	   
		for seqID in seqIDList[i]:
			seq=seqList[i][j]
			currentclass.append({'id': seqID, 'seq': seq})
			j=j+1
			count=count+1            
		taxadict[classname]=currentclass
		i=i+1
	#write to file
	with open(jsonfilename,"w") as json_file:
		json.dump(taxadict,json_file,encoding='latin1')

if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]
	#load data
	newfastafilename,seqids,sequences,taxa,classes,level = load_data(fastafilename,classificationfilename,classificationlevel)
	#represent sequences as matrix of k-mer frequencies
	filename=GetBase(fastafilename)
	matrixfilename=filename + "." + str(k) + ".matrix"
	command=path + "fasta2matrix.py " +  str(k) + " " + newfastafilename + " " + matrixfilename
	os.system(command)
	X,Y,nb_classes,input_length,data_max,seqIDList,seqList = load_matrix(matrixfilename,seqids,sequences,taxa,classes)
	#train data
	trainids=np.array(range(0,len(X)),dtype=int)
	traindata=X[trainids]
	trainlabels=Y[trainids]
	#training
	model = create_model(nb_classes,input_length)
	traindata = traindata.reshape(traindata.shape + (1,))
	trainlabels_bin=np_utils.to_categorical(trainlabels, nb_classes)
	model.fit(traindata, trainlabels_bin, nb_epoch=100, batch_size=20, verbose = 0)
	#save model
	modelname=filename.replace(".","_") + "_cnn_classifier"
	if level !="":
		modelname=filename + "_" + level + "_cnn_classifier"
	basename=modelname
	if "/" in modelname:
		basename=modelname[modelname.rindex("/")+1:]
	if os.path.isdir(modelname) == False:
		os.system("mkdir " + modelname)
	#save model
	classifiername=modelname + "/" + basename + ".classifier"
	model.save(classifiername)
	#save seqids for each classification
	jsonfilename=modelname + "/" + basename + ".classes"
	SaveClasses(jsonfilename,classes,seqIDList,seqList)
	#save config	
	configfilename=modelname + "/" + basename + ".config"
	SaveConfig(configfilename,classifiername,fastafilename,jsonfilename,classificationfilename,classificationlevel,k,data_max)
	
	

