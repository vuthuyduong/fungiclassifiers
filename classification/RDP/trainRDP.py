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


fastafilename=sys.argv[1]
classificationfilename=sys.argv[2] # 
classificationposition=sys.argv[3] 
rdpclassifierpath = sys.argv[4] # the path to the rdp classifier file classifier.jar

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

def GetSeqIDs(seqrecords):
	seqids=[]
	for seqrecord in seqrecords:
		seqids.append(seqrecord.id)
	return seqids

##############################################################################
# MAIN
##############################################################################
basefilename=GetBase(fastafilename)
#load seq records
seqrecords = list(SeqIO.parse(fastafilename, "fasta"))

#train the rdp classifier
rdpfastafilename= basefilename + ".rdp.fasta"
rdptaxaidfilename = basefilename +  ".rdp.tid"

#train the data using rdp model
trainseqids = GetSeqIDs(seqrecords)

#load taxonomic classifcation of the train dataset
species,genera,families,orders,classes,phyla,kingdoms,classifications,trainlabels=LoadClassification(trainseqids,classificationfilename,int(classificationposition))

#generate the taxaid file for the rdp model
GenerateTaxaIDs(species,genera,families,orders,classes,phyla,kingdoms,rdptaxaidfilename)

#generate fasta file for the rdp model
GenerateRDFFastaFile(trainseqids,classifications,fastafilename,rdpfastafilename)

trainfolder = basefilename + "_rdp_classifier"
if os.path.isdir(trainfolder) == False:
	os.system("mkdir " + trainfolder)
os.system("cp " + rdpclassifierpath + "/classifier/samplefiles/rRNAClassifier.properties " + trainfolder + "/")
traincommand="java -Xmx1g -jar " + rdpclassifierpath + "/classifier.jar train -o " + trainfolder + " -s " + rdpfastafilename + " -t " + rdptaxaidfilename
os.system(traincommand)





