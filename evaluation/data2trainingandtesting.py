#!/usr/bin/env python
from sklearn.model_selection import StratifiedKFold
import sys
import os
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Convolution1D
from keras.datasets import mnist
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution1D, MaxPooling1D
from keras.utils import np_utils
from keras import backend as K
import numpy as np

#parameters: sys.argv[1] = input dataset as matrix of k-mers
nome_train=sys.argv[1].split(".")[0]				
fold_number = 10
if len(sys.argv) >2:
	fold_number =int(sys.argv[2])

def load_data(file):
	classification=[]
	records= list(open(file, "r"))
	records=records[1:]
	X=[]
	for seq in records:
		elements=seq.split(",")
		X.append(elements[1:-1])
		level=elements[-1].split("\n")
		taxonname=level[0]
		classification.append(taxonname)
	taxaset=set(classification)
	classes=list(taxaset)
	
	Y=[]
	for taxonname in classification:
		Y.append(classes.index(taxonname))
	X=np.array(X,dtype=float)
	Y=np.array(Y,dtype=int)
	data_max= np.amax(X)
	X = X/data_max
	return X,Y,len(classes),len(X[0])

if __name__ == "__main__":
	n_folds = fold_number
	X,Y,nb_classes,input_length = load_data(sys.argv[1])
	kfold = StratifiedKFold(n_splits=n_folds, shuffle=True)
	if os.path.isdir("results") == False:
		os.system("mkdir results")
	i=1
	for train,test in kfold.split(X, Y):
                np.save("results/train_"+nome_train + "_"+str(i),train)
                np.save("results/test_"+nome_train + "_"+str(i),test)
                np.save("results/train_"+nome_train + "_"+str(i)+".labels",Y[train])
                np.save("results/test_"+nome_train + "_"+str(i)+".labels",Y[test])
		i=i+1
