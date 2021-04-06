#!/usr/bin/env python
# FILE: classifyRDP.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
if sys.version_info[0] >= 3:
	unicode = str
import math
import numpy as np
import os, argparse
from Bio import SeqIO
#from keras.utils import np_utils
import json
import random

parser=argparse.ArgumentParser(prog='classifyRDP.py',  
							   usage="%(prog)s [options] -i fastafile -c classifiername -cn classificationfilename -p classificationposition -rdp rdpclassifierpath -mp minproba -mc mincoverage -j variationjsonfilename",
							   description='''Script that classifies the sequences of the fasta files using RDP model. The classified sequences with a probability less than minproba will be verified using BLAST. The mincoverage is given for BLAST comparision. The json file is the the variation of the sequences within each group of the training dataset. This file is created during the training of the model, and used optionally for the verification of the classification. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-rdp','--rdpclassifierpath', required=True, help='the path of the RDP classifier classifier.jar.')
parser.add_argument('-o','--out', help='The folder name containing the model and associated files.') #optional
parser.add_argument('-c','--classifier', required=True, help='the folder containing the classifier.')

args=parser.parse_args()

testfastafilename= args.input
rdpclassifierpath = args.rdpclassifierpath
modelname=args.classifier
rdp_output=args.out


def GetBase(filename):
	if "." not in filename:
		return filename
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
	configfile=open(configfilename)
	for line in configfile:
		texts=line.split(": ")
		if texts[0]=="Classifier name":
			classifiername=texts[1].rstrip()
		if texts[0]=="Classification filename":
			classificationfilename=texts[1].rstrip()
		if texts[0]=="Column number to be classified":
			classificationpos=int(texts[1].rstrip())
		if texts[0]=="RDP fasta filename":
			rdpfastafilename=texts[1].rstrip()
		if texts[0]=="RDP ID filename":
			rdpidfilename=texts[1].rstrip()
		if texts[0]=="Classes filename":
			jsonfilename=texts[1].rstrip()	
	return jsonfilename,classifiername,classificationfilename,classificationpos,rdpfastafilename,rdpidfilename

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

def GetSeqIDs(seqrecords):
	seqids=[]
	for seqrecord in seqrecords:
		seqids.append(seqrecord.id)
	return seqids

def GetPredictedLabels(testseqids,rdpclassifiedfilename):
	predlabels =[""] * len(testseqids)
	scores=[0] * len(testseqids)
	ranks=[""] * len(testseqids)
	rdpfile=open(rdpclassifiedfilename)
	for line in  rdpfile:  
		words=line.split("\t")
		n=len(words)
		seqid = words[0]
		i=testseqids.index(seqid)
		scores[i]=float(words[n-1].rstrip())
		ranks[i]=words[n-2]
		predlabels[i]=words[n-3]
	rdpfile.close()
	return predlabels,scores,ranks

def SavePrediction(classdict,testseqids,testlabels,predlabels,probas,outputname):
	output=open(outputname,"w")
	output.write("SequenceID\tGiven classification\tPrediction\tFull classification of prediction\tProbability\n")
	i=0
	keys=classdict.keys()
	for i in range(0,len(testseqids)):
		seqid=testseqids[i]
		giventaxonname= testlabels[i]
		predictedname=predlabels[i]
		proba=probas[i]
		#predictedname=predictedname.encode('ascii', 'ignore')	
		classification=classdict[predictedname]['classification']
		output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + str(proba) + "\n")
		i=i+1
	output.close()
		
##############################################################################
# MAIN
##############################################################################
basefilename=GetBase(testfastafilename)


#load config of the model
classesfilename,classifiername,classificationfilename,classificationpos,rdpfastafilename,rdpidfilename=LoadConfig(modelname)

#load seq records
seqrecords = list(SeqIO.parse(testfastafilename, "fasta"))
testseqids = GetSeqIDs(seqrecords)

#load given taxonomic classifcation of the test dataset
testspecies,testgenera,testfamilies,testorders,testclasses,testphyla,testkingdoms,testclassifications,testlabels,rank=LoadClassification(testseqids,classificationfilename,classificationpos)

# run the model
basename=modelname
if "/" in basename:
	basename=basename[basename.rindex("/")+1:]

rdpclassifiedfilename = basefilename +  "." + basename + ".out"
if rdp_output==None or rdp_output=="":
	rdp_output=basefilename +  "." + basename + ".classified"


testcommand="java -Xmx1g -jar " + rdpclassifierpath + "/classifier.jar classify -t " + modelname + "/rRNAClassifier.properties -o " + rdpclassifiedfilename + " " + testfastafilename
os.system(testcommand)
testseqids = GetSeqIDs(seqrecords)
predlabels,probas,ranks=GetPredictedLabels(testseqids,rdpclassifiedfilename)
#load ref class dict
classdict={}
with open(classesfilename) as classesfile:
	classdict = json.load(classesfile)
SavePrediction(classdict,testseqids,testlabels,predlabels,probas,rdp_output)
print("The result is saved in the file: " + rdp_output)




