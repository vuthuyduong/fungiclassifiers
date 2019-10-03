## Training and classifying using the CNN model

To train the CNN model at the genus level, use the following command:

CNN/trainCNN.py traindataset.fas traindataset.classification 5 

where 5 is the position of the genera in the classification traindataset.classification file. Here the k-mer size is set to 6 as default. 

If we want to use another k-mer such as 8 for example, then use the following command:

CNN/trainCNN.py traindataset.fas traindataset.classification 5 8

To classify the test dataset using the classifier obtained after the training:

CNN/classifyCNN.py traindataset_genus_cnn_classifier testdataset.fas 

The result will be saved in the file testdataset.traindataset_genus_cnn_classifier.classified.

If we wish to calculate the BLAST similarity score for the predictions with a probability >= proba_threshold, then use the following command:

CNN/classifyCNN.py traindataset_genus_cnn_classifier testdataset.fas proba_threshold

## Training and classifying using the DBN model

To train the DBN model at the genus level, use the following command:

DBN/trainDBN.py traindataset.fas traindataset.classification 5 

where 5 is the position of the genera in the classification traindataset.classification file. Here the k-mer size is set to 6 as default. 

If we want to use another k-mer such as 8 for example, then use the following command:

DBN/trainDBN.py traindataset.fas traindataset.classification 5 8

To classify the test dataset using the classifier obtained after the training:

DBN/classifyDBN.py traindataset_genus_dbn_classifier testdataset.fas 

The result will be saved in the file testdataset.traindataset_genus_dbn_classifier.classified.

If we wish to calculate the BLAST similarity score for the predictions with a probability >= proba_threshold, then use the following command:

DBN/classifyDBN.py traindataset_genus_dbn_classifier testdataset.fas proba_threshold


## Training and classifying using the BLAST classification

To predict an optimal threshold for classifying the sequences at the genus level, for example, use the following command:

BLAST/trainBLAST.py traindataset.fas traindataset.classification 5 0.9 1 0.001 

where [0.9,1] is the interval of the optimal threshold, and 0.001 is the step used for searching. This command is used only when we don't know which BLAST threshold is used for classifying the sequences at the current taxonomic level.

To classify the sequences using the BLAST classification, used the following command:

BLAST/classifyBLAST.py testdataset.fas traindataset.fas traindataset.classification 5 0.97

where 0.97 is the optimal threshold that we use for classifying the sequences at the current level.

## Training and classifying using the RDP model

To train the RDP model at the genus level, use the following command:

RDP/trainRDP.py traindataset.fas traindataset.classification 5 path_to_the_file_classifer.jar

To classify the test dataset using the classifier obtained by the RDP training, use the following command:

RDP/classifyRDP.py RDP_classifier testdataset.fas path_to_the_file_classifer.jar

