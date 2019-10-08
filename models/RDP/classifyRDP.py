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
import json
import random


modelname = sys.argv[1]
rdpclassifierpath = sys.argv[2] # the path to the rdp classifier file classifier.jar
testfastafilename=sys.argv[3]
minprobaforBlast=1.1 #search for best match of the sequences by BLast only when minprobabilityforBlast <=1.0
if len(sys.argv)>4:
	minprobaforBlast=float(sys.argv[4])
classificationfilename=""
if len(sys.argv) > 5:
	classificationfilename=sys.argv[5] #
classificationposition=0
if len(sys.argv) >6: 
	classificationposition=sys.argv[6] 

mincoverage=300 #for ITS sequences, used only for comparing the sequences with BLAST

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
	data_max=0
	k=6
	configfile=open(configfilename)
	for line in configfile:
		texts=line.split(": ")
		if texts[0]=="Classifier name":
			classifiername=texts[1].rstrip()
		if texts[0]=="Classification filename":
			classificationfilename=texts[1].rstrip()
		if texts[0]=="Column number to be classified":
			classificationpos=int(texts[1].rstrip())
		if texts[0]=="Classes filename":
			jsonfilename=texts[1].rstrip()
	return jsonfilename,classifiername,classificationfilename,classificationpos

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
					orderclass=i
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


def SavePrediction(minprobaforBlast,mincoverage,classeswithsequences,testseqids,testlabels,predlabels,probas,outputname):
	output=open(outputname,"w")
	if minprobaforBlast <=1.0:
		output.write("Index\tSequenceID\tClassification\tPrediction\tProbability\tBLAST score\tBLAST sim\tBLAST coverage\tBLAST local sim\tBLAST local coverage\tNumber of sequences to be compared\tBest match ID\n")
	else:
		output.write("Index\tSequenceID\tClassification\tPrediction\tProbability\n")	
	for i in range(0,len(testseqids)):
		seqid=testseqids[i]
		giventaxonname= testlabels[i]
		predictedname=predlabels[i]
		proba=probas[i]
		if minprobaforBlast <= 1.0 and proba >= minprobaforBlast:
			bestrefseqID,score,sim,coverage,localsim,localcoverage,numberofsequences=ComputeBLASTScore(testseqrecords[i],predictedname,classeswithsequences,mincoverage)
			output.write(str(i) + "\t" + seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t" + str(proba) + "\t" + str(score) + "\t" + str(sim) + "\t" + str(coverage) + "\t" + str(localsim) + "\t" + str(localcoverage) + "\t" + str(numberofsequences) + "\t" + bestrefseqID + "\n")			
		else:
			output.write(str(i) + "\t" + seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t" + str(proba) + "\n")
		i=i+1
	output.close()
##############################################################################
# MAIN
##############################################################################
basefilename=GetBase(testfastafilename)


#load config of the model
jsonfilename,classifiername,classificationfilename,classificationposition=LoadConfig(modelname)

#load seq records
seqrecords = list(SeqIO.parse(testfastafilename, "fasta"))
testseqids = GetSeqIDs(seqrecords)

#load given taxonomic classifcation of the test dataset
testspecies,testgenera,testfamilies,testorders,testclasses,testphyla,testkingdoms,testclassifications,testlabels,rank=LoadClassification(testseqids,classificationfilename,int(classificationposition))

# run the model
basename=modelname
if "/" in basename:
	basename=basename[basename.rindex("/")+1:]

rdpclassifiedfilename = basefilename +  "." + basename + ".out"
rdp_output=basefilename +  "." + basename + ".classified"


testcommand="java -Xmx1g -jar " + rdpclassifierpath + "/classifier.jar classify -t " + modelname + "/rRNAClassifier.properties -o " + rdpclassifiedfilename + " " + testfastafilename
os.system(testcommand)
testseqids = GetSeqIDs(seqrecords)
predlabels,probas,ranks=GetPredictedLabels(testseqids,rdpclassifiedfilename)
classeswithsequences={}
testseqrecords=[]
if minprobaforBlast <=1.0:
	#load classes
	with open(jsonfilename) as json_file:
		classeswithsequences = json.load(json_file)
		testseqrecords=list(SeqIO.parse(testfastafilename, "fasta"))
print("The result is saved in the file: " + rdp_output)
SavePrediction(minprobaforBlast,mincoverage,classeswithsequences,testseqids,testlabels,predlabels,probas,rdp_output)





