#!/usr/bin/env python
# FILE: trainDBN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
import numpy as np
np.random.seed(1337)  # for reproducibility
#from dbn.tensorflow import SupervisedDBNClassification
import json
from Bio import SeqIO

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from sklearn.preprocessing import OneHotEncoder

parser=argparse.ArgumentParser(prog='trainDBN.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile -p classificationposition",
							   description='''Script that trains a DBN model to classify sequences''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The folder name containing the model and associated files.') #optional
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-k','--kmer', type=int, default=6, help='the k-mer for the representation of the sequences.')
parser.add_argument('-rank','--rank',default="species", help='the rank to load taxonomic names for the prediction.')


args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
rank=args.rank
k = args.kmer
modelname=args.out

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2] #taxonomy file 
#classificationlevel=int(sys.argv[3]) #the level of classification to get taxa from the taxonomy file
#k = 6
#if len(sys.argv) >4:
#	k= int(sys.argv[4])

class DBN:
	def __init__(self, hidden_layers_structure, learning_rate_rbm=0.05, learning_rate=0.1,
                 n_epochs_rbm=10, n_iter_backprop=500, batch_size=32, activation_function='relu',
                 dropout_p=0.1, verbose=False):
		self.hidden_layers_structure = hidden_layers_structure
		self.learning_rate_rbm = learning_rate_rbm
		self.learning_rate = learning_rate
		self.n_epochs_rbm = n_epochs_rbm
		self.n_iter_backprop = n_iter_backprop
		self.batch_size = batch_size
		self.activation_function = activation_function
		self.dropout_p = dropout_p
		self.verbose = verbose

	def fit(self, X_train, y_train):
		# Preprocess target labels
		enc = OneHotEncoder()
		y_train_encoded = enc.fit_transform(y_train.reshape(-1, 1)).toarray()

		# Initialize DBN model
		model = Sequential()

		# Add hidden layers
		for units in self.hidden_layers_structure:
			model.add(Dense(units, activation=self.activation_function))
			model.add(Dropout(self.dropout_p))

		# Add output layer
		model.add(Dense(y_train_encoded.shape[1], activation='softmax'))

		# Compile model
		model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=self.learning_rate),
		              loss='categorical_crossentropy',
		              metrics=['accuracy'])

		# Train the model
		model.fit(X_train, y_train_encoded, epochs=self.n_iter_backprop, batch_size=self.batch_size,
		          verbose=self.verbose)

		self.model = model

	def predict(self, X_test):
		# Make predictions
		return self.model.predict(X_test)

	def predict_classes(self, X_test):
		# Make class predictions
		return tf.argmax(self.model.predict(X_test), axis=1)

	def save_model(self, filepath):
		self.model.save(filepath)
        
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

def GetClassificationPos(classificationfilename,rank):
	classificationfile=open(classificationfilename, errors='ignore')
	header=next(classificationfile)
	texts=header.split("\t")
	i=0
	classificationpos=-1
	for text in texts:
		text=text.rstrip()
		if text.lower()==rank.lower():
			classificationpos=i
		i=i+1	
	classificationfile.close()
	if classificationpos <0:
		print("Please specify a rank for building the model.")
		os.system(exit)	
	return classificationpos
	
def load_data(fastafilename,classificationfilename,rank):
	#load seqrecords
	seqids=[]
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load classification
	classificationfile= open(classificationfilename, errors='ignore')
	header=next(classificationfile)
	texts=header.split("\t")
	classnames=[]
	newseqrecords=[]
	classdict={}
	header=""
	classificationpos=GetClassificationPos(classificationfilename,rank)
	for line in classificationfile:
		texts=line.split("\t")
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
				if len(classification) > len(classdict[classname]['classification']):
					classdict[classname]['classification']=classification
				classdict[classname]['seqids'].append(seqid)
			else:	
				classdict.setdefault(classname,{})
				classdict[classname]['classification']=classification
				classdict[classname]['seqids']=[seqid]
	classificationfile.close()
	newfastafilename=GetBase(fastafilename) + "." + str(classificationpos) + ".fasta"
	#write selected sequences to file
	SeqIO.write(newseqrecords, newfastafilename, "fasta")
	return newfastafilename,seqids,classnames,classdict

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
		index=list(classdict.keys()).index(taxonname)
		Y.append(index)
		X.append(elements[1:])
	X=np.array(X,dtype=float)
	Y=np.array(Y,dtype=int)
	#normalize
	data_max= np.amax(X)
	X = X/data_max
	return X,Y,len(set(classnames)),len(X[0]),data_max

# def create_model():
#  	classifier = SupervisedDBNClassification(hidden_layers_structure=[256, 256],
#                                               learning_rate_rbm=0.05,
#                                               learning_rate=0.1,
#                                               n_epochs_rbm=10,
#                                               n_iter_backprop=500,
#                                               batch_size=32,
#                                               activation_function='relu',
#                                               dropout_p=0.1,verbose=False)
#  	return classifier

# def create_model(input_length):
#     model = Sequential([
#         Dense(256, activation='relu', input_shape=(input_length,1)),  # Assuming input_shape is defined
#         Dense(256, activation='relu'),
#         Dropout(0.1),  # Add dropout layer with 0.1 dropout rate
#         Dense(nb_classes, activation='softmax')  # Assuming nb_classes is defined
#     ])
#     model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.1),
#                   loss='categorical_crossentropy',
#                   metrics=['accuracy'])
#     return model

def create_model():
 	classifier = DBN(hidden_layers_structure=[256, 256],
                                              learning_rate_rbm=0.05,
                                              learning_rate=0.1,
                                              n_epochs_rbm=10,
                                              n_iter_backprop=500,
                                              batch_size=32,
                                              activation_function='relu',
                                              dropout_p=0.1,verbose=False)
 	return classifier

def SaveConfig(configfilename,classifiername,fastafilename,jsonfilename,classificationfilename,rank,kmer,data_max):
	if not classifiername.startswith("/"):
		classifiername=os.getcwd() + "/" + classifiername
	if not fastafilename.startswith("/"):
		fastafilename=os.getcwd() + "/" + fastafilename
	if not jsonfilename.startswith("/"):
		jsonfilename=os.getcwd() + "/" + jsonfilename
	if not classificationfilename.startswith("/"):
		classificationfilename=os.getcwd() + "/" + classificationfilename
	model="dbn"
	#save the config:classifierfilename, model, classificationfilename,classificationpos,k-mer
	configfile=open(configfilename,"w")
	configfile.write("Classifier name: " + classifiername + "\n")
	configfile.write("Model: " + model + "\n")
	configfile.write("Data max: " + str(data_max) + "\n")
	configfile.write("K-mer number: " + str(kmer) + "\n")
	configfile.write("Fasta filename: " + fastafilename + "\n")
	configfile.write("Classification filename: " + classificationfilename + "\n")
	configfile.write("Rank to be classified: " + rank + "\n")
	configfile.write("Classes filename: " + jsonfilename + "\n")
	configfile.close()	

# def GetLevel(classificationfilename,classificationpos):
# 	classificationfile=open(classificationfilename)
# 	header=next(classificationfile)
# 	if header.startswith("#"):
# 		level=header.split("\t")[classificationpos]				     
# 		level=level.rstrip()
# 	else:
# 		level=str(classificationpos)   
# 	return level

def SaveClasses(jsonfilename,classdict):
	#write json dict
	with open(jsonfilename,"w") as json_file:
		json.dump(classdict,json_file)
		
if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]
	filename=GetBase(fastafilename)
	if modelname==None or modelname=="":
		modelname=filename.replace(".","_") + "_dbn_classifier" 
		if rank !="":
			modelname=filename + "_" + rank + "_dbn_classifier"
	basename=modelname
	if "/" in modelname:
		basename=modelname[modelname.rindex("/")+1:]
	#load data
	newfastafilename,seqids,classnames,classdict = load_data(fastafilename,classificationfilename,rank)
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
	model = None # Clearing the NN.
	model = create_model()
	model.fit(traindata, trainlabels)
	#save model
	#modelname=filename.replace(".","_") + "_dbn_classifier"
	
	if os.path.isdir(modelname) == False:
		os.system("mkdir " + modelname)
	#save model
	classifiername=modelname + "/" + basename + ".classifier"
	model.save_model(classifiername)
	#save seqids for each classification
	jsonfilename=modelname + "/" + basename + ".classes"
	SaveClasses(jsonfilename,classdict)
	#save classification file:
	os.system("cp " + classificationfilename + " " + modelname)
	#save config	
	configfilename=modelname + "/" + basename + ".config"
	SaveConfig(configfilename,classifiername,fastafilename,jsonfilename,classificationfilename,rank,k,data_max)
	print("The classifier is saved in the folder " + modelname + ".")
	

