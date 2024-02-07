## Training and classifying using the CNN model

To train the CNN model at the genus level, use the following command:

*CNN/trainCNN.py -i traindataset.fas -c traindataset.classification -rank species* 

where the classification file will have a format like the .classifiction file in the data folder. Here the k-mer size is set to 6 as default. 

If we want to use another k-mer such as 8 for example, then use the following command:

*CNN/trainCNN.py -i traindataset.fas -c traindataset.classification -rank species -k 6*

To classify the test dataset using the classifier obtained after the training:

*CNN/classifyCNN.py -c traindataset_genus_cnn_classifier -i testdataset.fas*

The result will be saved in the file testdataset.traindataset_genus_cnn_classifier.out.

## Training and classifying using the DBN model

To train the DBN model at the genus level, use the following command:

*DBN/trainDBN.py -i traindataset.fas -c traindataset.classification -rank species* 

Here the k-mer size is set to 6 as default. 

If we want to use another k-mer such as 8 for example, then use the following command:

*DBN/trainDBN.py -i traindataset.fas -c traindataset.classification -rank species 5 -k 6*

To classify the test dataset using the classifier obtained after the training:

*DBN/classifyDBN.py -c traindataset_genus_dbn_classifier -i testdataset.fas*

The result will be saved in the file testdataset.traindataset_genus_dbn_classifier.out.

## Training and classifying using the RDP model

To train the RDP model at the genus level, use the following command:

*RDP/trainRDP.py -i traindataset.fas -c traindataset.classification -p 5 -rdp path_to_the_file_classifer.jar*

To classify the test dataset using the classifier obtained by the RDP training, use the following command:

*RDP/classifyRDP.py -c RDP_classifier -i testdataset.fas -rdp path_to_the_file_classifer.jar*


