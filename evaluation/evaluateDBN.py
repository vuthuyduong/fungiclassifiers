#!/usr/bin/env python
# FILE: evaluateDBN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import os
import numpy as np
np.random.seed(1337)  # for reproducibility
from sklearn.datasets import load_digits
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics.classification import accuracy_score
from dbn.tensorflow import SupervisedDBNClassification

from datetime import datetime, date
from datetime import timedelta

matrixfilename=sys.argv[1] # the matrix obtained by fasta2matrix combined with taxonomic classification
datapath=sys.argv[2] #deeplearing results'path

def create_model():
	classifier = SupervisedDBNClassification(hidden_layers_structure=[256, 256],
                                             learning_rate_rbm=0.05,
                                             learning_rate=0.1,
                                             n_epochs_rbm=10,
                                             n_iter_backprop=500,
                                             batch_size=32,
                                             activation_function='relu',
                                             dropout_p=0.1,verbose=False)
	return classifier

def train_and_evaluate_model (model, X_train, Y_train, X_test, Y_test):
	beginning=datetime.now()
	model.fit(X_train, Y_train)
	end=datetime.now()
	t1=(end-beginning).total_seconds()
	beginning=datetime.now()
	Y_pred = model.predict(X_test)
	classifying_score = accuracy_score(Y_test, Y_pred)
	end=datetime.now()
	t2=(end-beginning).total_seconds()
    	return Y_pred, Y_test,t1,t2




def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(matrixfilename):
	records= list(open(matrixfilename, "r"))
	records=records[1:]
	classification=[]
	seqnames=[]
	length=0
	X=[]
	for seq in records:
		elements=seq.split(",")
		X.append(elements[1:-1])
		level=elements[-1].split("\n")
		taxonname=level[0]
		classification.append(taxonname)
		seqnames.append(elements[0])
	taxaset=set(classification)
	classes=list(taxaset)
	Y=[]
	for taxonname in classification:
		Y.append(classes.index(taxonname))
	X=np.array(X,dtype=float)
	Y=np.array(Y,dtype=int)
	data_max= np.amax(X)
	X = X/data_max
	return X,Y,classification,seqnames,len(classes),len(X[0])



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

reportfilename=GetBase(matrixfilename) + ".dbn.report"


#Load classes, classification:
X,Y,classification,seqnames,nb_classes,input_length= LoadClassification(matrixfilename)


#generate fasta files for training and test dataset
filenames = os.listdir(datapath)
prefix = GetBase(matrixfilename) 
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
	if (filename.startswith(train_prefix)) and ("labels" not in filename):
		trainids = np.load(datapath+"/"+filename)
	if (filename.startswith(test_prefix)) and ("labels" in filename):
		if not basefilename.replace("test_","") in datasets:
			datasets.append(basefilename.replace("test_","").replace(".labels",""))			


#generate report
report=open(reportfilename,"w")
report.write("Test dataset\tNumber of sequences in the test dataset\tTrain dataset\tNumber of sequences in the train dataset\tAccuracy of DBN\tNumber of unidentified sequences by DBN\tTraining time\tClassifying time\n")
for dataset in datasets:
	print(dataset)
	testdataset=datapath+"/test_"+dataset
	traindataset=datapath+"/train_"+dataset
	#load testids, trainids:
	testids=np.load(datapath+"/test_"+dataset+".npy")
	trainids=np.load(datapath+"/train_"+dataset+".npy")

	#load labels 
	train_labels=np.load(datapath +"/train_" + dataset +".labels.npy")
	test_labels=np.load(datapath +"/test_" + dataset +".labels.npy")
	

	#train the data
	model = None # Clearing the NN.
	model = create_model()
	pred, Y_test,t1,t2 = train_and_evaluate_model(model, X[trainids], Y[trainids], X[testids], Y[testids])
	np.save(datapath + "/preds_" + dataset + ".labels",pred)

	#load the predicted labels of the test dataset
	preds_labels=np.load(datapath+"/preds_"+dataset+".labels.npy")

	#get all labels
	alllabels=GetAllLabels(trainids,train_labels,testids,test_labels)

	#calulate accuracy for deeplearning
	acc_of_deeplearning,unidentifiedlist_by_deeplearning=ComputeAccuracy(train_labels,testids,test_labels,preds_labels)

	#save prediction by deeplearning
	deeplearning_output=datapath+"/test_"+dataset+".dbn.out"
	SavePrediction(classification,seqnames,alllabels,testids,preds_labels,deeplearning_output)
	#Print the result
	report.write(testdataset + "\t" + str(len(testids)) + "\t" + traindataset + "\t" + str(len(trainids)) + "\t" + str(acc_of_deeplearning)+ "\t" + str(len(unidentifiedlist_by_deeplearning)) + "\t" + str(t1) + "\t" + str(t2) + "\n")
	print(testdataset + "\t" + str(len(testids)) +  "\t" + traindataset + "\t" + str(len(trainids)) + "\t" + str(acc_of_deeplearning)+ "\t" + str(len(unidentifiedlist_by_deeplearning)) + "\t" + str(t1) + "\t" + str(t2) +"\n")

report.close()
