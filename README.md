# Fungiclassifiers

Fungiclassifiers is the official implementation of the classifiers using two different models namely convolutional neural network (CNN) (LeCun et al. 2015) and deep belief network (DBN) (Hinton & Salakhutdinov 2006, Hinton et al. 2006) for fungal classification.

## Implementation

Fungiclassifiers consists of four parts: Classification, Evaluation, Data and classifiers.

The training, classifying, and evaluation of the BLAST, RDP, CNN and DBN models were implemented in Python 2.7. We used the Keras library (www.keras.io) with tensorflow backend for CNN and the code available at https://github.com/albertbup/deep-belief-network for DBN, as used in (Fiannaca et al. 2018, https://github.com/IcarPA-TBlab/MetagenomicDC). The RDP Clas-sifier was downloaded from https://github.com/rdpstaff/classifier.  For the BLAST classification, BLAST version 2.6.0 was installed. The fMLC tool to cluster the sequences for predicting optimal thresholds for taxonomic classification were downloaded from https://github.com/FastMLC/fMLC. Source code and datasets are available at https://github.com/vuthuyduong/fungiclassifiers.

## Classification

[Classification](https://github.com/vuthuyduong/fungiclassifiers/tree/master/classification) consists of the scripts for training and  for training and classifying a data

## Evaluation

[evaluation](https://github.com/vuthuyduong/fungiclassifiers/tree/master/evaluation)

## Data
There are two datasets available as inputs for fungiclassifiers. The "small" dataset contains ~4000 ITS yeast sequences, checked and validated by the specialists at the Westerdijk Fungal Biodiversity Institute. This dataset were analyzed and released in [Vu D. et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5192050/). The "large" dataset contains ~350K ITS fungal sequences downloaded from GenBank (https://www.ncbi.nlm.nih.gov/) which was used in [Vu D. et al. 2014](https://www.nature.com/articles/srep06837) to evaluate the speed of MLC.

[yeast dataset](http://www.westerdijkinstitute.nl/Download/SmallDatasetOf4KYeastITSSequences.zip) the small demo dataset.

[top50 UNITE](https://github.com/vuthuyduong/fungiclassifiers/tree/master/data/top50UNITE) the large demo dataset. 

## Classifiers
[Classifiers](https://github.com/vuthuyduong/fungiclassifiers/tree/master/classifiers)



## Contact person 

Duong Vu (d.vu@wi.knaw.nl)

## References
LeCun, Y., Bengio, Y., Hinton, G. (2015). Deep learning. Nature, 521, 436-44.
Hinton, G.E., Salakhutdinov, R.R (2006). Reducing the Dimensionality of Data with Neural Networks. Science 313, 504–7.



