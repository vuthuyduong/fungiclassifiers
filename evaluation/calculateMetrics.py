#!/usr/bin/env python
# FILE: calculateMetrics.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import numpy as np
import os


from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
#from keras.layers import Dense
#from keras.layers import Convolution1D

datapath=sys.argv[1]
method=sys.argv[2] #cnn,dbn,rdp,blast

def GetBase(filename):
	if "." in filename:
		return filename[:-(len(filename)-filename.rindex("."))]
	else:
		return filename

def CalculateMetrics(test_labels,pred_labels,labels,allabels): 
	accuracy=accuracy_score(test_labels,pred_labels)
	precision,recall,fscore,support=precision_recall_fscore_support(test_labels,pred_labels,average='macro')
	precisionvector,recallvector,fscorevector,support=precision_recall_fscore_support(test_labels,pred_labels,labels=labels)
	mcc=matthews_corrcoef(test_labels,pred_labels)
	cohenkappascore=cohen_kappa_score(test_labels,pred_labels,labels=labels)
	confusionmatrix=confusion_matrix(test_labels,pred_labels,labels=alllabels)
	return accuracy,precision,recall,fscore,precisionvector,recallvector,fscorevector,cohenkappascore,mcc,confusionmatrix

def GetLabels(classifiedfilename):
	classifiedfile=open(classifiedfilename)
	next(classifiedfile)
	seqids=[]
	seqnames=[]
	test_labels=[]
	pred_labels=[]
	for line in classifiedfile:
		words=line.rstrip().split("\t")
		seqids.append(words[0])
		seqnames.append(words[1])
		testlabel=""
		if len(words) >2:
			testlabel=words[2]
		predlabel=""
		if len(words) >3:
			predlabel=words[3]
		proba=1
#		if len(words) >4:
#			proba=float(words[4])
#			if proba < 0.8:
#				predlabel=""
		test_labels.append(testlabel)
		pred_labels.append(predlabel)
	classifiedfile.close()	
	return seqids,seqnames,test_labels,pred_labels

def Filter(train_labels,test_labels,testlabels,predlabels):
	filteredtestlabels=[]
	filteredpredlabels=[]
	i=0
	for test_label in test_labels:
		if test_label in train_labels:
			filteredtestlabels.append(testlabels[i])
			filteredpredlabels.append(predlabels[i])
		i=i+1
	return filteredtestlabels,filteredpredlabels

##############################################################################
# MAIN
##############################################################################

reportfilename=GetBase(datapath) + "." + method +  ".metrics"
matrixfilename=GetBase(datapath) + "." + method +  ".matrix"

#generate fasta files for training and test dataset
filenames = os.listdir(datapath)
train_prefix = "train_" 
test_prefix="test_"
pred_predix="preds_"

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
report.write("Test dataset\tNumber of sequences in the test dataset\tTrain dataset\tNumber of sequences in the train dataset\tAccuracy\tRecall\tPrecision\tf_score\tCohen Kappa score\tMcc\tfiltered accuracy\tfiltered recall\tfiltered precision\tfiltered fscore\tcohen kappa score\tfilteredMcc\tNumber of unidentified sequences\tNumber of labels in the train dataset\n")
precisions=[]
recalls=[]
fscores=[]
testnumbers=[]
alltestlabels=[]
numberoftrainsequences=[]

filteredprecisions=[]
filteredrecalls=[]
filteredfscores=[]
filteredtestnumbers=[]
allfilteredtestlabels=[]
filterednumberoftrainsequences=[]
confusionmatrics=[]
alllabels=[]
for dataset in datasets:
	testdataset=datapath+"/test_"+dataset
	classifiedfilename=testdataset + "." + method + ".out"
	if os.path.isfile(classifiedfilename)==False:
		classifiedfilename=testdataset + "." + method + ".classified"
	if os.path.isfile(classifiedfilename)==False and method.lower()=="rdp":
		classifiedfilename=testdataset + "_by_rdp.out"	
	if os.path.isfile(classifiedfilename)==False:
		continue
	#get test and pred labels
	seqids,seqnames,testlabels,predlabels= GetLabels(classifiedfilename)
	labels=list(set(testlabels+predlabels))
	alllabels=alllabels+labels
alllabels=list(set(alllabels))
alllabels.sort()	
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

	unidentifiedlist=[]
	numberoftrainlabels=len(set(train_labels))
	
	for i in range(0,len(testids)):
		if test_labels[i] not in train_labels:
			unidentifiedlist.append(testids[i])
		
	#load the predicted labels of the test dataset
	#preds_labels=np.load(datapath+"/preds_"+ dataset + "." + method + ".labels.npy")

	#get all labels
	#alllabels=GetAllLabels(trainids,train_labels,testids,test_labels)
	classifiedfilename=testdataset + "." + method + ".out"
	if os.path.isfile(classifiedfilename)==False:
		classifiedfilename=testdataset + "." + method + ".classified"
	if os.path.isfile(classifiedfilename)==False and method.lower()=="rdp":
		classifiedfilename=testdataset + "_by_rdp.out"	
	if os.path.isfile(classifiedfilename)==False:
		continue
	#get test and pred labels
	seqids,seqnames,testlabels,predlabels= GetLabels(classifiedfilename)
	labels=list(set(testlabels+predlabels))
	
	#filter the test labels that have a label in the train dataset
	filteredtestlabels,filteredpredlabels=Filter(train_labels,test_labels,testlabels,predlabels)
	filteredlabels=list(set(filteredtestlabels + filteredpredlabels))
	
	#calculate number of sequences of the (filtered) labels in the train dataset
	numberoftrainsequencesvector=[]
	for label in labels:
		trainsequencenumber=0
		if label in testlabels:
			labelid=test_labels[testlabels.index(label)]
			trainsequencenumber=list(train_labels).count(labelid)
		numberoftrainsequencesvector.append(trainsequencenumber)
	filterednumberoftrainsequencesvector=[]
	for label in filteredlabels:
		trainsequencenumber=0
		if label in testlabels:
			labelid=test_labels[testlabels.index(label)]
			trainsequencenumber=list(train_labels).count(labelid)
		filterednumberoftrainsequencesvector.append(trainsequencenumber)
	
	
	#calulate metrics for all the test labels
	labels.sort()
	accuracy,precision,recall,fscore,precisionvector,recallvector,fscorevector,cohenkappascore,mcc,confusionmatrix=CalculateMetrics(testlabels,predlabels,labels,alllabels)
	#calulate metrics for all the test labels that have a label in the train dataset
	filteredlabels.sort()
	filteredaccuracy,filteredprecision,filteredrecall,filteredfscore,filteredprecisionvector,filteredrecallvector,filteredfscorevector,filteredcohenkappascore,filteredmcc,filteredconfusionmatrix=CalculateMetrics(filteredtestlabels,filteredpredlabels,filteredlabels,alllabels)
	i=0
	for label in labels:
		if label in alltestlabels:
			index=alltestlabels.index(label)
			testnumbers[index]=testnumbers[index]+1
			precisions[index]=precisions[index] + precisionvector[i]
			recalls[index]=recalls[index] + recallvector[i]
			fscores[index]=fscores[index] + fscorevector[i]
			numberoftrainsequences[index]=numberoftrainsequences[index] + numberoftrainsequencesvector[i]
			
		else:
			alltestlabels.append(label)
			testnumbers.append(1)
			precisions.append(precisionvector[i])
			recalls.append(recallvector[i])
			fscores.append(fscorevector[i])
			numberoftrainsequences.append(numberoftrainsequencesvector[i])
		i=i+1	
	i=0
	for label in filteredlabels:
		if label in allfilteredtestlabels:
			index=allfilteredtestlabels.index(label)
			filteredtestnumbers[index]=filteredtestnumbers[index]+1
			filteredprecisions[index]=filteredprecisions[index] + filteredprecisionvector[i]
			filteredrecalls[index]=filteredrecalls[index] + filteredrecallvector[i]
			filteredfscores[index]=filteredfscores[index] + filteredfscorevector[i]
			filterednumberoftrainsequences[index]=filterednumberoftrainsequences[index] + filterednumberoftrainsequencesvector[i]
		else:
			allfilteredtestlabels.append(label)
			filteredtestnumbers.append(1)
			filteredprecisions.append(filteredprecisionvector[i])
			filteredrecalls.append(filteredrecallvector[i])
			filteredfscores.append(filteredfscorevector[i])
			filterednumberoftrainsequences.append(filterednumberoftrainsequencesvector[i])
		i=i+1		
	if confusionmatrics==[]:
		confusionmatrics=filteredconfusionmatrix
	else:
		confusionmatrics=confusionmatrics+filteredconfusionmatrix
	
	#Print the result
	report.write(testdataset + "\t" + str(len(testids)) + "\t" + traindataset + "\t" + str(len(trainids)) + "\t" + str(accuracy) + "\t" + str(recall)  + "\t" + str(precision) + "\t" + str(fscore)  + "\t" + str(cohenkappascore)  + "\t" + str(mcc) + "\t" + str(filteredaccuracy) + "\t" + str(filteredrecall)+ "\t" + str(filteredprecision) + "\t" + str(filteredfscore)  + "\t" + str(filteredcohenkappascore) + "\t" + str(filteredmcc) + "\t" + str(len(unidentifiedlist)) + "\t" + str(numberoftrainlabels) + "\n")
	print(testdataset + "\t" + str(len(testids)) +  "\t" + traindataset + "\t" + str(len(trainids)) + "\t" +str(accuracy) + "\t" + str(recall)+ "\t" + "\t" + str(precision) + "\t" + str(fscore)  + "\t" + str(cohenkappascore) + "\t" + str(mcc) + "\t" + str(filteredaccuracy) + "\t" + str(filteredrecall) + "\t" + str(filteredprecision) + "\t" + str(filteredfscore)  +  "\t" + str(filteredcohenkappascore) + "\t" + str(filteredmcc) + "\t" + str(len(unidentifiedlist)) + "\t" + str(numberoftrainlabels) + "\n")

report.write("\n")
report.write("Average metrics for each class that has a label in the train dataset: " + "\n")
report.write("Class name\tRecall\tPrecision\tfscore\tNumber of train sequences\n")
i=0
for label in allfilteredtestlabels:
	j=allfilteredtestlabels.index(label)
	n=filteredtestnumbers[i]
	precision=float(filteredprecisions[i])/n
	recall=float(filteredrecalls[i])/n
	fscore=float(filteredfscores[i])/n
	seqnumber=int(round(float(filterednumberoftrainsequences[i])/n,0))
	report.write(label + "\t" + str(recall) + "\t" + str(precision) + "\t" + str(fscore) + "\t" + str(seqnumber) + "\n" )
	i=i+1
	
	
print("The accuracy of the model is saved in the file " + reportfilename )
report.close()

#generate matrix
#print(filteredconfusionmatrix)
matrix=open(matrixfilename,"w")
for label in alllabels:
	matrix.write("\t" + label)
matrix.write("\n")	
for i in range(0,len(alllabels)):
	matrix.write(alllabels[i])
	for j in range(0,len(alllabels)):
		matrix.write("\t" + str(confusionmatrics[i][j]))
	matrix.write("\n")	
matrix.close()
print("The confusion matrix of the model is saved in the file " + matrixfilename )
