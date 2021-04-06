#!/usr/bin/env python
# FILE: trainRDP.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
if sys.version_info[0] >= 3:
	unicode = str
#import math
#import numpy as np
import os, argparse
from Bio import SeqIO
#from keras.utils import np_utils
import json
#import random

parser=argparse.ArgumentParser(prog='trainRDP.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile,-p classificationposition -rdp rdpclassifierpath",
							   description='''Script that trains a RDP model to classify sequences''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The folder name containing the model and associated files.') #optional
parser.add_argument('-c','--classification', help='the classification file in tab. format.')
parser.add_argument('-p','--classificationpos', type=int,required=True,default=0, help='the classification position to load the classification.')
parser.add_argument('-rdp','--rdpclassifierpath', required=True, help='the path of the RDP classifier classifier.jar.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
classificationpos=args.classificationpos
rdpclassifierpath = args.rdpclassifierpath
modelname=args.out

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2] # 
#classificationposition=sys.argv[3] 
#rdpclassifierpath = sys.argv[4] # the path to the rdp classifier file classifier.jar

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]


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
	kingdompos=1
	phylumpos=2 #the position of the phyla in the classification file
	classpos=3
	orderpos=4
	familypos=5
	genuspos=6
	speciespos=7
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
			if pos==0:
				pos=speciespos
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

def GenerateRDFFastaFile(seqids,labels,classifications,trainfastafilename,rdpfastafilename):
	classdict={}
	rdpfastafile = open(rdpfastafilename,"w")
	trainfastafile =open(trainfastafilename)
	seqid=""
	seq=""
	for line in trainfastafile:
		if line.startswith(">"):
			#add the previous sequences to classes
			seqid=line.rstrip().replace(">","")
			i=seqids.index(seqid)
			label=labels[i]
			classification=classifications[i]
			if label in classdict.keys():
				if len(classification) > len(classdict[label]['classification']):
					classdict[label]['classification']=classification
				classdict[label]['seqids'].append(seqid)
			else:
				classdict.setdefault(label,{})
				classdict[label]['classification']=classification
				classdict[label]['seqids']=[seqid]
			rdpfastafile.write(">" + seqid + "\t" + classification + "\n")
		else:
			seq=seq + line.rstrip()
			rdpfastafile.write(line)	
	rdpfastafile.close()
	trainfastafile.close()
	return classdict
	
def GetSeqIDs(seqrecords):
	seqids=[]
	for seqrecord in seqrecords:
		seqids.append(seqrecord.id)
	return seqids

def SaveClasses(jsonfilename,classdict):
	#write json dict
	with open(jsonfilename,"w") as json_file:
		json.dump(classdict,json_file,encoding='latin1')

def SaveConfig(configfilename,classifiername,fastafilename,jsonfilename,classificationfilename,classificationpos,rdpfastafilename,rdpidfilename):
	if not classifiername.startswith("/"):
		classifiername=os.getcwd() + "/" + classifiername
	if not fastafilename.startswith("/"):
		fastafilename=os.getcwd() + "/" + fastafilename
	if not classificationfilename.startswith("/"):
		classificationfilename=os.getcwd() + "/" + classificationfilename
	model="rdp"
	#save the config:classifierfilename, model, classificationfilename,classificationpos,k-mer
	configfile=open(configfilename,"w")
	configfile.write("Model: " + model + "\n")
	configfile.write("Fasta filename: " + fastafilename + "\n")
	configfile.write("Classification filename: " + classificationfilename + "\n")
	configfile.write("Column number to be classified: " + str(classificationpos) + "\n")
	configfile.write("Classes filename: " + jsonfilename + "\n")
	configfile.write("RDP fasta filename: " + rdpfastafilename + "\n")
	configfile.write("RDP ID filename: " + rdpidfilename + "\n")
	configfile.close()

##############################################################################
# MAIN
##############################################################################
basefilename=GetBase(fastafilename)
#load seq records
seqrecords = list(SeqIO.parse(fastafilename, "fasta"))
#train the data using rdp model
trainseqids = GetSeqIDs(seqrecords)

#load taxonomic classifcation of the train dataset
species,genera,families,orders,classes,phyla,kingdoms,classifications,trainlabels,level=LoadClassification(trainseqids,classificationfilename,classificationpos)

#save model
if modelname==None or modelname=="":
	modelname = basefilename + "_rdp_classifier"
	if level !="":
		modelname=basefilename + "_" + level + "_rdp_classifier"
if os.path.isdir(modelname) == False:
	os.system("mkdir " + modelname)	
basename=modelname
if "/" in modelname:
	basename=modelname[modelname.rindex("/")+1:]	
#basename=basefilename + "_" + level + "_rdp_classifier"
if "/" in basename:
	basename = basename[len(basename)-basename.rindex("/"):]
os.system("cp " + rdpclassifierpath + "/classifier/samplefiles/rRNAClassifier.properties " + modelname + "/")

#train the rdp classifier
rdpfastafilename= modelname + "/" + basefilename + ".rdp.fasta"
rdptaxaidfilename = modelname + "/" + basefilename +  ".rdp.tid"

#generate the taxaid file for the rdp model
GenerateTaxaIDs(species,genera,families,orders,classes,phyla,kingdoms,rdptaxaidfilename)

#generate fasta file for the rdp model
classdict=GenerateRDFFastaFile(trainseqids,trainlabels,classifications,fastafilename,rdpfastafilename)


traincommand="java -Xmx1g -jar " + rdpclassifierpath + "/classifier.jar train -o " + modelname + " -s " + rdpfastafilename + " -t " + rdptaxaidfilename
os.system(traincommand)
#save seqids for each classification
jsonfilename=modelname + "/" + basename + ".classes"
SaveClasses(jsonfilename,classdict)
#save config	
configfilename=modelname + "/" + basename + ".config"
SaveConfig(configfilename,modelname,fastafilename,jsonfilename,classificationfilename,classificationpos,rdpfastafilename,rdptaxaidfilename)
print("The classifier is saved in the folder " + modelname + ".")
	
	




