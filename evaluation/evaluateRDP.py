#!/usr/bin/env python
# FILE: trainingresults2report.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import math
import numpy as np
import os
from Bio import SeqIO
#from keras.utils import np_utils
from datetime import datetime, date
from datetime import timedelta


fastafilename=sys.argv[1]
classificationfilename=sys.argv[2] # 
classificationposition=sys.argv[3] 
datapath=sys.argv[4] #the folder containing the train and test datasets in .npy format
rdpclassifierpath = sys.argv[5] # the path to the rdp classifier file classifier.jar

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(seqids,classificationfilename,pos):
	kingdompos=0
	phylumpos=1 #the position of the phyla in the classification file
	classpos=2
	orderpos=3
	familypos=4
	genuspos=5
	speciespos=6
	species=[""] * len(seqids)
	genera=[""] * len(seqids)
	families=[""] * len(seqids)
	orders=[""] * len(seqids)
	classes=[""] * len(seqids)
	phyla=[""] * len(seqids)
	kingdoms=[""] * len(seqids)
	classificationfile=open(classificationfilename)
	classifications=[""] * len(seqids)
	labels=[""] * len(seqids)
	level= 0 #kingdom
	sep="\t"
	for line in classificationfile:
		if ";" in line:
			sep=";"
		words=line.split(sep)
		if line.startswith("#"):			
			i=0
			for word in words:
				word=word.rstrip()
				if word=="species":
					speciespos=i
				if word=="genus":
					genuspos=i
				if word=="family":
					familypos=i
				if word=="order":
					orderpos=i
				if word=="class":
					orderclass=i
				if word=="phylum":
					phylumpos=i
				if word=="kingdom":
					kingdompos=i
			if pos==speciespos:
				level=6
			if pos==genuspos:
				level=5
			if pos==familypos:
				level=4
			if pos==orderpos:
				level=3
			if pos==classpos:
				level=2
			if pos==phylumpos:
				level=1	
			if pos==kingdompos:
				level=1			
			continue 
		seqid=words[0].replace(">","").rstrip()
		if seqid in seqids:
			index=seqids.index(seqid)
			classification="Root"
			kingdom="Fungi"
			if level >=0:
				if kingdompos >0:
					kingdom=words[kingdompos].rstrip()
				kingdoms[index]= kingdom
				classification=classification + ";" + kingdom
			if level >=1:
				phyla[index] = words[phylumpos].rstrip()
				classification=classification + ";" + words[phylumpos].rstrip()
			if level >=2:
				classes[index]= words[classpos].rstrip()
				classification=classification + ";" +  words[classpos].rstrip()
			if level >=3:
				orders[index]=words[orderpos].rstrip()
				classification=classification + ";" + words[orderpos].rstrip()
			if level>=4:
				families[index]=words[familypos].rstrip()
				classification=classification + ";" + words[familypos].rstrip()
			if level>=5:
				genera[index]=words[genuspos].rstrip()
				classification=classification + ";" + words[genuspos].rstrip()
			if level>=6:
				species[index]=words[speciespos].rstrip()
				classification=classification + ";" + words[speciespos].rstrip()
			classifications[index]=classification
			labels[index]=words[pos].rstrip()
	return species,genera,families,orders,classes,phyla,kingdoms,classifications,labels

def GenerateTaxaIDs(species,genera,families,orders,classes,phyla,kingdoms,taxaidfilename):
	classificationfile=open(classificationfilename)
	taxids=[]
	taxa=[]
	parenttaxids=[]
	depths=[]
	ranks=[]
	#add root
	taxids.append(0)
	taxa.append("Root")
	parenttaxids.append(-1)
	depths.append(0)
	ranks.append("rootrank")
	i=0
	for kingdom in kingdoms:
		#add  kingdom
		kingdom_id=len(taxa)
		if kingdom !="":
			if kingdom in taxa:
				kingdom_id=taxa.index(kingdom)
			else:
				taxids.append(len(taxa))
				taxa.append(kingdom)
				parenttaxids.append(0)
				depths.append(1)
				ranks.append("domain")
		phylum_id=len(taxa)
		if len(phyla)>0:
			if phyla[i] in taxa:
				phylum_id=taxa.index(phyla[i])
			elif phyla[i]!= "":
				taxids.append(len(taxa))
				taxa.append(phyla[i])
				parenttaxids.append(kingdom_id)
				depths.append(2)
				ranks.append("phylum")
		class_id=len(taxa)
		if len(classes)>i:
			if classes[i] in taxa:
				class_id=taxa.index(classes[i])
			elif classes[i] != "":
				taxids.append(len(taxa))
				taxa.append(classes[i])
				parenttaxids.append(phylum_id)
				depths.append(3)
				ranks.append("class")
		order_id=len(taxa)
		if len(orders)>i:
			if orders[i] in taxa:
				order_id=taxa.index(orders[i])
			elif orders[i] != "":
				taxids.append(len(taxa))
				taxa.append(orders[i])
				parenttaxids.append(class_id)
				depths.append(4)
				ranks.append("order")
		family_id=len(taxa)
		if len(families) >i:
			if families[i] in taxa:
				family_id=taxa.index(families[i])
			elif families[i] != "":
				taxids.append(len(taxa))
				taxa.append(families[i])
				parenttaxids.append(order_id)
				depths.append(5)
				ranks.append("family")
		genus_id=len(taxa)
		if len(genera) > i:
			if genera[i] in taxa:
				genus_id=taxa.index(genera[i])
			elif genera[i] != "":
				taxids.append(len(taxa))
				taxa.append(genera[i])
				parenttaxids.append(family_id)
				depths.append(6)
				ranks.append("genus")
		species_id=len(taxa)
		if len(species) > i:
			if species[i] in taxa:
				species_id=taxa.index(species[i])
			elif species[i] != "":
				taxids.append(len(taxa))
				taxa.append(species[i])
				parenttaxids.append(genus_id)
				depths.append(7)
				ranks.append("species")
		i=i+1
	#write to taxaidfilename
	taxaidfile=open(taxaidfilename,"w")
	i=0
	for taxid in taxids:
		taxaidfile.write(str(taxid) + "*" + taxa[i] + "*" + str(parenttaxids[i]) + "*" + str(depths[i]) + "*" +ranks[i]+"\n")
		i=i+1
	taxaidfile.close()

def GenerateRDFFastaFile(seqids,classifications,trainfastafilename,rdpfastafilename):
	rdpfastafile = open(rdpfastafilename,"w")
	trainfastafile =open(trainfastafilename)
	for line in trainfastafile:
		if line.startswith(">"):
			seqid=line.rstrip().replace(">","")
			rdpfastafile.write(">" + seqid + "\t" + classifications[seqids.index(seqid)] + "\n")
		else:
			rdpfastafile.write(line)
	rdpfastafile.close()
	trainfastafile.close()

def GenerateFastafile(prefixfilename,ids,seqrecords):	
	seqfile = open(prefixfilename + ".fasta","w")
	i=0
	for seqrecord in seqrecords:
		if i in ids:
			#save this sequence in the seqfile
			seqfile.write(">" + seqrecord.description + "\n")
			seqfile.write(str(seqrecord.seq) + "\n")
		i=i+1	
	seqfile.close()

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

def GetSeqIDs(ids,seqrecords):
	seqids=[]
	for id in ids:
		seqids.append(seqrecords[id].id)
	return seqids


def GetPredictedLabels(testseqids,rdpclassifiedfilename):
	predlabels =[""] * len(testseqids)
	scores=[0] * len(testseqids)
	ranks=[""] * len(testseqids)
	rdpfile=open(rdpclassifiedfilename)
	for line in  rdpfile:  
		line=line.split("\n")[0]
		words=line.split("\t")
		n=len(words)
		seqid = words[0]
		i=testseqids.index(seqid)
		scores[i]=float(words[n-1])
		ranks[i]=words[n-2]
		predlabels[i]=words[n-3]
	rdpfile.close()
	return predlabels,scores,ranks

def SavePrediction(testids,testseqids,testlabels,predlabels,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tClassification\tPrediction\tn")
	i=0
	for id in testids:
		output.write(str(id) + "\t" + testseqids[i] + "\t" + testlabels[i] + "\t"  + predlabels[i] + "\n")
		i=i+1
	output.close()
##############################################################################
# MAIN
##############################################################################

reportfilename=GetBase(fastafilename) + ".rdp.report"

path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]
#load seq records
seqrecords = list(SeqIO.parse(fastafilename, "fasta"))

#generate fasta files for training and test dataset
filenames = os.listdir(datapath)
prefix = GetBase(fastafilename) 
if "/" in prefix:
	prefix =prefix[prefix.rindex("/")+1:]
train_prefix = "train_" + prefix
test_prefix="test_" + prefix
pred_predix="preds_" + prefix

datasets=[]

for filename in filenames:
	if not ".npy" in filename:
		continue 
	basefilename = GetBase(filename)
	if (filename.startswith(test_prefix)) and ("tanh" not in filename):
		testids = np.load(datapath+"/"+filename)
		GenerateFastafile(datapath+"/"+basefilename,testids,seqrecords)
		if not basefilename.replace("test_","") in datasets:
			datasets.append(basefilename.replace("test_",""))
	if (filename.startswith(train_prefix)) and ("tanh" not in filename):
		trainids = np.load(datapath+"/"+filename)
		GenerateFastafile(datapath+"/"+basefilename,trainids,seqrecords)		
#generate report
report=open(reportfilename,"w")
report.write("Test dataset\tNumber of sequences in the test dataset\tTrain dataset\tNumber of sequences in the train dataset\tAccuracy of RDP model\tNumber of unidentified sequences by RDP model\tTraining time\tClassifying time\n")
for dataset in datasets:
	print(dataset)
	#load testids, trainids:
	testids=np.load(datapath+"/test_"+dataset+".npy")
	trainids=np.load(datapath+"/train_"+dataset+".npy")
	testdataset=datapath+"/test_"+dataset+".fasta"
	traindataset=datapath+"/train_"+dataset+".fasta"
	rdpfastafilename= datapath+"/train_"+dataset + ".rdp.fasta"
	rdptaxaidfilename = datapath+"/train_"+dataset + ".rdp.tid"
	rdpclassifiedfilename = datapath+"/test_"+dataset + ".rdp.txt"

	#train the data using rdp model
	trainseqids = GetSeqIDs(trainids,seqrecords)

	#load taxonomic classifcation of the train dataset
	species,genera,families,orders,classes,phyla,kingdoms,classifications,trainlabels=LoadClassification(trainseqids,classificationfilename,int(classificationposition))
		
	#generate the taxaid file for the rdp model
	GenerateTaxaIDs(species,genera,families,orders,classes,phyla,kingdoms,rdptaxaidfilename)

	#generate fasta file for the rdp model
	GenerateRDFFastaFile(trainseqids,classifications,traindataset,rdpfastafilename)

	#train the rdp classifier
	beginning=datetime.now()
	trainfolder = datapath + "/train_"+ dataset + "_rdp_classifier" 
	if os.path.isdir(trainfolder) == False:
		os.system("mkdir " + trainfolder)
	os.system("cp " + rdpclassifierpath + "/classifier/samplefiles/rRNAClassifier.properties " + trainfolder + "/")
	traincommand="java -Xmx1g -jar " + rdpclassifierpath + "/classifier.jar train -o " + trainfolder + " -s " + rdpfastafilename + " -t " + rdptaxaidfilename
	os.system(traincommand)
	end=datetime.now()
	t1=(end-beginning).total_seconds()
	#classify the test dataset
	beginning=datetime.now()
	testcommand="java -Xmx1g -jar " + rdpclassifierpath + "/classifier.jar classify -t " + trainfolder + "/rRNAClassifier.properties -o " + rdpclassifiedfilename + " " + testdataset
	os.system(testcommand)
	end=datetime.now()
	t2=(end-beginning).total_seconds()
	#load classification result
	testseqids = GetSeqIDs(testids,seqrecords)
	predlabels,scores,ranks=GetPredictedLabels(testseqids,rdpclassifiedfilename)	

	#load given taxonomic classifcation of the test dataset
	testspecies,testgenera,testfamilies,testorders,testclasses,testphyla,testkingdoms,testclassifications,testlabels=LoadClassification(testseqids,classificationfilename,int(classificationposition))

	#Caculate the accuracy based on rdp
	acc_of_rdp,unidentifiedlist_by_rdp=ComputeAccuracy(trainlabels,testids,testlabels,predlabels)
	report.write(testdataset + "\t" + str(len(testids)) + "\t" + traindataset + "\t" + str(len(trainids)) + "\t" + str(acc_of_rdp)+ "\t" + str(len(unidentifiedlist_by_rdp)) + "\t" + str(t1) + "\t" + str(t2) + "\n")
	#Save prediction by rdp
	rdp_output=datapath+"/test_"+dataset+"_by_rdp.out"
	SavePrediction(testids,testseqids,testlabels,predlabels,rdp_output)
	print(testdataset + "\t" + str(len(testids)) +  "\t" + traindataset + "\t" + str(len(trainids)) + "\t" + str(acc_of_rdp)+ "\t" + str(len(unidentifiedlist_by_rdp)) + "\t" + str(t1) + "\t" + str(t2) + "\n")

report.close()
