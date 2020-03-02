#!/usr/bin/env python
# FILE: evaluateDBN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import os
import numpy as np
np.random.seed(1337)  # for reproducibility
from sklearn.metrics import precision_recall_fscore_support
from sklearn.datasets import load_digits
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics.classification import accuracy_score
from dbn.tensorflow import SupervisedDBNClassification

from datetime import datetime, date
from datetime import timedelta

matrixfilename=sys.argv[1] # the matrix obtained by fasta2matrix combined with taxonomic classification
datapath=sys.argv[2] #deeplearing results'path
reportfilename=""
if len(sys.argv) > 3:
	reportfilename=sys.argv[3]
evalnumber=0 #evaluate all datasets at once if evalnumber=0, individually otherwise
if len(sys.argv) >4:
	evalnumber=int(sys.argv[4])

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
	#classifying_score = accuracy_score(Y_test, Y_pred)
	end=datetime.now()
	t2=(end-beginning).total_seconds()
	probas = model.predict_proba(X_test)
	return Y_pred, probas, Y_test,t1,t2

def GetBase(filename):
	base=filename
	if "." in filename:
		base=filename[:-(len(filename)-filename.rindex("."))]
	if "/" in base:
		base=base[base.rindex("/")+1:]
	return base 
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

def CalculateMetrics(train_labels,testids,test_labels,pred_labels): 
	precision,recall,fscore,support=precision_recall_fscore_support(test_labels,pred_labels,average='micro')
	return precision,recall,fscore

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

def SavePrediction(classification,seqnames,alllabels,testids,pred_labels,probas,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tClassification\tPrediction\tProbability\n")
	i=0
	for id in testids:
		proba =max(probas[i])
		predictedname=""
		if pred_labels[i] in alllabels:
			j=alllabels.index(pred_labels[i])
			predictedname=classification[j]
		output.write(str(id) + "\t" + str(seqnames[id]) + "\t" + classification[id] + "\t"  + predictedname + "\t" + str(proba) +"\n")
		i=i+1
	output.close()
##############################################################################
# MAIN
##############################################################################
if reportfilename=="":
	reportfilename=GetBase(matrixfilename) + ".dbn.report"


#Load classes, classification:
X,Y,classification,seqnames,nb_classes,input_length= LoadClassification(matrixfilename)


#generate fasta files for training and test dataset
filenames = os.listdir(datapath)
prefix = GetBase(matrixfilename) 
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
	if (filename.startswith(train_prefix)) and ("labels" not in filename):
		trainids = np.load(datapath+"/"+filename)
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
		dataset=testdataset.split("test_")[1]
		print(dataset)
		evaluateddatasets.append(dataset)
	report.close()
	reportexists=True
#generate report
if reportexists==False:
	report=open(reportfilename,"w")
	report.write("Test dataset\tNumber of sequences in the test dataset\tTrain dataset\tNumber of sequences in the train dataset\tNumber of unidentified sequences\tTraining time\tClassifying time\tCreation time\n")
	report.close()
#write to report
i=0
for dataset in datasets:
	if dataset in evaluateddatasets:
		continue
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
	beginning=datetime.now()
	model = None # Clearing the NN.
	model = create_model()
	end=datetime.now()
	creationtime=(end-beginning).total_seconds()
	pred, probas, Y_test,t1,t2 = train_and_evaluate_model(model, X[trainids], Y[trainids], X[testids], Y[testids])
	np.save(datapath + "/preds_" + dataset + ".dbn.labels",pred)

	#load the predicted labels of the test dataset
	preds_labels=np.load(datapath+"/preds_"+ dataset + ".dbn.labels.npy")

	#get all labels
	alllabels=GetAllLabels(trainids,train_labels,testids,test_labels)
	
	#calulate metrics
	#precision,recall,fscore=CalculateMetrics(train_labels,testids,test_labels,preds_labels)
	
	#get unidentified list
	unidentifiedlist=[]
	for i in range(0,len(testids)):
		if test_labels[i] not in train_labels:
			unidentifiedlist.append(testids[i])
			
	#save prediction by deeplearning
	deeplearning_output=datapath+"/test_"+dataset+".dbn.classified"
	SavePrediction(classification,seqnames,alllabels,testids,preds_labels,probas,deeplearning_output)
	#Print the result
	report=open(reportfilename,"a")
	report.write(testdataset + "\t" + str(len(testids)) + "\t" + traindataset + "\t" + str(len(trainids))  + "\t" + str(len(unidentifiedlist)) + "\t" + str(t1) + "\t" + str(t2) + "\t" + str(creationtime) + "\n")
	report.close()
	print(testdataset + "\t" + str(len(testids)) +  "\t" + traindataset + "\t" + str(len(trainids))  + "\t" + str(len(unidentifiedlist)) + "\t" + str(t1) + "\t" + str(t2) + str(creationtime) +"\n")
	i=i+1
	if i==evalnumber:
		break
print("The running time of the model is saved in the file " + reportfilename + ". The results of the classification are saved as .dbn.classified files in the folder " + datapath + "." )
#report.close()
