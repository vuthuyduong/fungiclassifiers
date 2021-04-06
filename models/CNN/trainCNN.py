#!/usr/bin/env python
# FILE: trainCNN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Convolution1D, MaxPooling1D
from keras.layers import Dropout, Activation, Flatten
from keras.utils import np_utils
import numpy as np
import json
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='trainCNN.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile -p classificationposition",
							   description='''Script that trains a CNN model to classify sequences''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The folder name containing the model and associated files.') #optional
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-p','--classificationpos', required=True, type=int, default=0, help='the classification position to load the classification.')
parser.add_argument('-k','--kmer', type=int, default=6, help='the k-mer for the representation of the sequences.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
classificationpos=args.classificationpos
k = args.kmer
modelname=args.out

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2] #taxonomy file 
#classificationlevel=int(sys.argv[3]) #the level of classification to get taxa from the taxonomy file
#k = 6
#if len(sys.argv) >4:
#	k= int(sys.argv[4])

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def GetTaxonomicClassification(classificationpos,header,texts):
	classification=""
	p_s=len(texts)
	p_g=len(texts)
	p_f=len(texts)
	p_o=len(texts)
	p_c=len(texts)
	p_p=len(texts)
	p_k=len(texts)
	i=0
	for text in header.split("\t"):
		text=text.rstrip()
		if text.lower()=="species":
			p_s=i
		elif text.lower()=="genus":
			p_g=i	
		elif text.lower()=="family":
			p_f=i	
		elif text.lower()=="order":
			p_o=i	
		elif text.lower()=="class":
			p_c=i	
		elif text.lower()=="phylum":
			p_p=i	
		elif text.lower()=="kingdom":
			p_k=i	
		i=i+1 
	species="s__"
	genus="g__"
	family="f__"
	order="o__" 
	bioclass="c__"
	phylum="p__"
	kingdom="k__"
	if p_s< len(texts):
		species="s__" + texts[p_s].rstrip()
	if p_g< len(texts):
		genus="g__" + texts[p_g].rstrip()	
	if p_f< len(texts):
		family="f__" + texts[p_f].rstrip()	
	if p_o< len(texts):
		order="o__" + texts[p_o].rstrip()
	if p_c< len(texts):
		bioclass="c__" + texts[p_c].rstrip()	
	if p_p< len(texts):
		phylum="p__" + texts[p_p].rstrip()
	if p_k< len(texts):
		kingdom="k__" + texts[p_k].rstrip()	
	if classificationpos==p_s:
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family + ";"+ genus+";"+species
	elif classificationpos==p_g:
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family + ";"+ genus
	elif classificationpos==p_f:
		classification=kingdom +";"+phylum + ";"+bioclass +";"+order+";"+family 
	elif classificationpos==p_o:
		classification=kingdom +";"+phylum + ";"+bioclass + ";"+order
	elif classificationpos==p_c:
		classification=kingdom +";"+phylum + ";"+bioclass 
	elif classificationpos==p_p:
		classification=kingdom +";"+phylum
	elif classificationpos==p_k:
		classification=kingdom
	else:
		classification=texts[classificationpos]
	return classification
		
def load_data(modelname,fastafilename,classificationfilename,classificationpos):
	#load seqrecords
	seqids=[]
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load classification
	classificationfile= open(classificationfilename)
	classnames=[]
	level=""
	newseqrecords=[]
	classdict={}
	header=""
	for line in classificationfile:
		texts=line.split("\t")
		if line.startswith("#"):
			header=line
			if classificationpos < len(texts):
				level=texts[classificationpos].rstrip()
			continue 		
		seqid=texts[0].replace(">","").rstrip()
		if not seqid in seqrecords.keys():
			continue
		classname=""
		if classificationpos < len(texts):
			classname=texts[classificationpos].rstrip()
			#classname=unicode(classname,errors='ignore')
		if classname !="":
			newseqrecords.append(seqrecords[seqid])
			seqids.append(seqid)
			classnames.append(classname)
			classification=GetTaxonomicClassification(classificationpos,header,texts)
			if classname in classdict.keys():
				if len(classification) > classdict[classname]['classification']:
					classdict[classname]['classification']=classification
				classdict[classname]['seqids'].append(seqid)
			else:	
				classdict.setdefault(classname,{})
				classdict[classname]['classification']=classification
				classdict[classname]['seqids']=[seqid]
	classificationfile.close()
	basename=GetBase(fastafilename)
	if "/" in basename:
		basename=basename[basename.rindex("/")+1:]
	newfastafilename=modelname + "/" + basename + "." + str(classificationpos) + ".fasta"
	#write selected sequences to file
	SeqIO.write(newseqrecords, newfastafilename, "fasta")
	return newfastafilename,seqids,classnames,classdict,level

def load_matrix(matrixfilename,seqids,classnames,classdict):
	#load matrix
	vectors= list(open(matrixfilename, "r"))
	vectors=vectors[1:]
	X=[]
	Y=[]
	for vector in vectors:
		elements=vector.split(",")
		seqid=elements[0]
		taxonname= classnames[seqids.index(seqid)]
		index=classdict.keys().index(taxonname)
		Y.append(index)
		X.append(elements[1:])
	X=np.array(X,dtype=float)
	Y=np.array(Y,dtype=int)
	data_max=0
	#data_max= np.amax(X)
	#X = X/data_max
	return X,Y,len(set(classnames)),len(X[0]),data_max

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

def GetLevel(classificationfilename,classificationpos):
	classificationfile=open(classificationfilename)
	header=next(classificationfile)
	if header.startswith("#"):
		level=header.split("\t")[classificationpos]				     
		level=level.rstrip()
	else:
		level=str(classificationpos)   
	return level

def SaveClasses(jsonfilename,classdict):
	#write json dict
	with open(jsonfilename,"w") as json_file:
		json.dump(classdict,json_file,encoding='latin1')

if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]
	filename=GetBase(fastafilename)
	level=GetLevel(classificationfilename,classificationpos)
	if modelname==None or modelname=="":
		modelname=filename.replace(".","_") + "_cnn_classifier" 
		if level !="":
			modelname=filename + "_" + level + "_cnn_classifier"
	basename=modelname
	if "/" in modelname:
		basename=modelname[modelname.rindex("/")+1:]
	#load data
	newfastafilename,seqids,classnames,classdict,level = load_data(modelname,fastafilename,classificationfilename,classificationpos)
	#represent sequences as matrix of k-mer frequencies
	matrixfilename=filename + "." + str(k) + ".matrix"
	command=path + "fasta2matrix.py " +  str(k) + " " + newfastafilename + " " + matrixfilename
	os.system(command)
	X,Y,nb_classes,input_length,data_max = load_matrix(matrixfilename,seqids,classnames,classdict)
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
#	modelname=filename.replace(".","_") + "_cnn_classifier"
	
	if os.path.isdir(modelname) == False:
		os.system("mkdir " + modelname)
	#save model
	classifiername=modelname + "/" + basename + ".classifier"
	model.save(classifiername)
	#save seqids for each classification
	jsonfilename=modelname + "/" + basename + ".classes"
	SaveClasses(jsonfilename,classdict)
	#save config	
	configfilename=modelname + "/" + basename + ".config"
	SaveConfig(configfilename,classifiername,fastafilename,jsonfilename,classificationfilename,classificationpos,k,data_max)
	print("The classifier is saved in the folder " + modelname + ".")
	

