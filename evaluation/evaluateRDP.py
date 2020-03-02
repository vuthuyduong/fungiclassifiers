#!/usr/bin/env python
# FILE: evaluateRDP.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import math
import numpy as np
import os
from Bio import SeqIO
#from keras.utils import np_utils
from sklearn.metrics import precision_recall_fscore_support
from datetime import datetime, date
from datetime import timedelta

fastafilename=sys.argv[1]
classificationfilename=sys.argv[2] # 
classificationposition=sys.argv[3] 
datapath=sys.argv[4] #the folder containing the train and test datasets in .npy format
rdpclassifierpath = sys.argv[5] # the path to the rdp classifier file classifier.jar
reportfilename=""
if len(sys.argv) > 6:
	reportfilename=sys.argv[6]

def GetBase(filename):
	base=filename
	if "." in filename:
		base=filename[:-(len(filename)-filename.rindex("."))]
	if "/" in base:
		base=base[base.rindex("/")+1:]
	return base 

def GetUnidentifiedSpeciesName(seqid):
	name="unidentified_species_of_" + seqid
	return name

def GetUnidentifiedGenusName(seqid,species):
	name=""
	if species!="":
		name="unidentified_genus_of_" + species
	else:
		name="unidentified_genus_of_" + seqid
	return name

def GetUnidentifiedFamilyName(seqid,species,genus):
	name=""
	if genus!="":
		name="unidentified_family_of_" + genus
	elif species!="":
		name="unidentified_family_of_" + species
	else:
		name="unidentified_family_of_" + seqid
	return name

def GetUnidentifiedOrderName(seqid,species,genus,family):
	name=""
	if family !="":
		name="unidentified_order_of_" + family
	elif genus!="":
		name="unidentified_order_of_" + genus
	elif species!="":
		name="unidentified_order_of_" + species
	else:
		name="unidentified_order_of_" + seqid
	return name

def GetUnidentifiedClassName(seqid,species,genus,family,order):
	name=""
	if order!="":
		name="unidentified_class_of_" + order
	elif family !="":
		name="unidentified_class_of_" + family
	elif genus!="":
		name="unidentified_class_of_" + genus
	elif species!="":
		name="unidentified_class_of_" + species
	else:
		name="unidentified_class_of_" + seqid
	return name

def GetUnidentifiedPhylumName(seqid,species,genus,family,order,bioclass):
	name=""
	if bioclass!="":
		name="unidentified_phylum_of_" + bioclass
	elif order!="":
		name="unidentified_phylum_of_" + order
	elif family !="":
		name="unidentified_phylum_of_" + family
	elif genus!="":
		name="unidentified_phylum_of_" + genus
	elif species!="":
		name="unidentified_phylum_of_" + species
	else:
		name="unidentified_phylum_of_" + seqid
	return name

def GetUnidentifiedKingdomName(seqid,species,genus,family,order,bioclass,phylum):
	name=""
	if phylum!="":
		name="unidentified_kingdom_of_" + phylum
	elif bioclass!="":
		name="unidentified_kingdom_of_" + bioclass
	elif order!="":
		name="unidentified_kingdom_of_" + order
	elif family !="":
		name="unidentified_kingdom_of_" + family
	elif genus!="":
		name="unidentified_kingdom_of_" + genus
	elif species!="":
		name="unidentified_kingdom_of_" + species
	else:
		name="unidentified_kingdom_of_" + seqid
	return name
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
	#load species, genera, families, orders, classes, phyla, kingdoms for sequences
	for line in classificationfile:
		line=unicode(line,errors='ignore')
		if ";" in line:
			sep=";"
		words=line.split(sep)
		if line.startswith("#"):			
			i=0
			for word in words:
				word=word.lower().rstrip()
				if word=="species":
					speciespos=i
				if word=="genus":
					genuspos=i
				if word=="family":
					familypos=i
				if word=="order":
					orderpos=i
				if word=="class":
					classpos=i
				if word=="phylum":
					phylumpos=i
				if word=="kingdom":
					kingdompos=i
				i=i+1
			continue 
		kingdom=""
		phylum=""
		bioclass=""
		order=""
		family=""
		genus=""
		spec=""
		if kingdompos>0 and len(words)>kingdompos:
			kingdom=words[kingdompos].rstrip()
		if phylumpos>0 and len(words)>phylumpos:
			phylum=words[phylumpos].rstrip()
		if classpos>0 and len(words)>classpos:
			bioclass=words[classpos].rstrip()
		if orderpos>0 and len(words)>orderpos:
			order=words[orderpos].rstrip()
		if familypos>0 and len(words)>familypos:
			family=words[familypos].rstrip()
		if genuspos>0 and len(words)>genuspos:
			genus=words[genuspos].rstrip()
		if speciespos>0 and len(words)>speciespos:
			spec=words[speciespos].rstrip()
		seqid=words[0].replace(">","").rstrip()
		if kingdom=="":
			kingdom=GetUnidentifiedKingdomName(seqid,spec,genus,family,order,bioclass,phylum)
		if phylum=="":
			phylum=GetUnidentifiedPhylumName(seqid,spec,genus,family,order,bioclass)
		if bioclass=="":
			bioclass=GetUnidentifiedClassName(seqid,spec,genus,family,order)
		if order=="":
			order=GetUnidentifiedOrderName(seqid,spec,genus,family)
		if family=="":
			family=GetUnidentifiedFamilyName(seqid,spec,genus)
		if genus=="":
			genus=GetUnidentifiedGenusName(seqid,spec)
		if spec=="":
			spec=GetUnidentifiedSpeciesName(seqid)
		label=""
		if pos==speciespos:
			level=6
			label=spec
			rank="species"
		if pos==genuspos:
			level=5
			label=genus
			rank="genus"
		if pos==familypos:
			level=4
			label=family
			rank="family"
		if pos==orderpos:
			level=3
			label=order
			rank="order"
		if pos==classpos:
			level=2
			label=bioclass
			rank="class"
		if pos==phylumpos:
			level=1	
			label=phylum
			rank="phylum"
		if pos==kingdompos:
			level=1	
			label=kingdom
			rank="kingdom"
		if seqid in seqids:
			index=seqids.index(seqid)
			classification="Root"
			kingdom="Fungi"
			if level >=0:
				kingdoms[index]= kingdom
				classification=classification + ";" + kingdom
			if level >=1:
				phyla[index] = phylum
				classification=classification + ";" + phylum
			if level >=2:
				classes[index]= bioclass
				classification=classification + ";" +  bioclass
			if level >=3:
				orders[index]=order
				classification=classification + ";" + order
			if level>=4:
				families[index]=family
				classification=classification + ";" + family
			if level>=5:
				genera[index]=genus
				classification=classification + ";" + genus
			if level>=6:
				species[index]=spec
				classification=classification + ";" + spec
			classifications[index]=classification
			labels[index]=label
	return species,genera,families,orders,classes,phyla,kingdoms,classifications,labels,rank

def GenerateTaxaIDs(species,genera,families,orders,classes,phyla,kingdoms,taxaidfilename):
	#classificationfile=open(classificationfilename)
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

def CalculateMetrics(train_labels,testids,test_labels,pred_labels): 
	precision,recall,fscore,support=precision_recall_fscore_support(test_labels,pred_labels,average='micro')
	return precision,recall,fscore

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

def SavePrediction(testids,testseqids,testlabels,predlabels,probas,outputname):
	output=open(outputname,"w")
	output.write("Sequence index\tSequenceID\tClassification\tPrediction\tProbability\n")
	i=0
	for id in testids:
		proba=probas[i]
		output.write(str(id) + "\t" + testseqids[i] + "\t" + testlabels[i] + "\t"  + predlabels[i] + "\t" + str(proba) + "\n")
		i=i+1
	output.close()
##############################################################################
# MAIN
##############################################################################
if reportfilename=="":
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
if "." in prefix:
	prefix=prefix[0:prefix.index(".")]	
train_prefix = "train_" + prefix
test_prefix="test_" + prefix
pred_predix="preds_" + prefix

datasets=[]
print(test_prefix)
for filename in filenames:
	if not ".npy" in filename:
		continue 
	basefilename = GetBase(filename)
	if (filename.startswith(test_prefix)) and ("labels" not in filename):
		testids = np.load(datapath+"/"+filename)
		GenerateFastafile(datapath+"/"+basefilename,testids,seqrecords)
		if not basefilename.replace("test_","") in datasets:
			datasets.append(basefilename.replace("test_",""))
	if (filename.startswith(train_prefix)) and ("labels" not in filename):
		trainids = np.load(datapath+"/"+filename)
		GenerateFastafile(datapath+"/"+basefilename,trainids,seqrecords)		

#generate report
report=open(reportfilename,"w")
report.write("Test dataset\tNumber of sequences in the test dataset\tTrain dataset\tNumber of sequences in the train dataset\tNumber of unidentified sequences\tTraining time\tClassifying time\n")
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
	species,genera,families,orders,classes,phyla,kingdoms,classifications,trainlabels,rank=LoadClassification(trainseqids,classificationfilename,int(classificationposition))
		
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
	predlabels,probas,ranks=GetPredictedLabels(testseqids,rdpclassifiedfilename)	

	#load given taxonomic classifcation of the test dataset
	testspecies,testgenera,testfamilies,testorders,testclasses,testphyla,testkingdoms,testclassifications,testlabels,rank=LoadClassification(testseqids,classificationfilename,int(classificationposition))

	#calulate metrics
	#precision,recall,fscore=CalculateMetrics(trainlabels,testids,testlabels,predlabels)
	
	#get unidentified list
	unidentifiedlist=[]
	for i in range(0,len(testids)):
		if testlabels[i] not in trainlabels:
			unidentifiedlist.append(testids[i])
			
	#Save prediction by rdp
	rdp_output=datapath+"/test_"+dataset+".rdp.classified"
	SavePrediction(testids,testseqids,testlabels,predlabels,probas,rdp_output)
	report.write(testdataset + "\t" + str(len(testids)) + "\t" + traindataset + "\t" + str(len(trainids))  + "\t" + str(len(unidentifiedlist)) + "\t" + str(t1) + "\t" + str(t2) + "\n")
	print(testdataset + "\t" + str(len(testids)) +  "\t" + traindataset + "\t" + str(len(trainids)) + "\t" + str(len(unidentifiedlist)) + "\t" + str(t1) + "\t" + str(t2) +"\n")
print("The running time of the model is saved in the file " + reportfilename + ". The results of the classification are saved as .rdp.classified files in the folder " + datapath + "." )
report.close()
