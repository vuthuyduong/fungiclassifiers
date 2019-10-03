## Training and classifying using the CNN model

To train the CNN model at the genus level, use the following command:

CNN/trainCNN.py traindataset.fas traindataset.classification 5 

where 6 is the position of the genera in the classification traindataset.classification file. Here the k-mer size is set to 6 as default. 

If we want to use another k-mer such as 8 for example, then use the following command:

CNN/trainCNN.py traindataset.fas traindataset.classification 6 8

To classify the test dataset using the classifier obtained after the training:

CNN/classifyCNN.py
