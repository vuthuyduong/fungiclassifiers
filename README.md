# Fungiclassifiers

Fungiclassifiers is the official implementation of the classifiers two different models namely convolutional neural network 
(CNN) (LeCun et al. 2015) and deep belief network (DBN) (Hinton & Salakhutdinov 2006, Hinton et al. 2006) for fungal classification.

## Implementation



## Classification

[classification](https://github.com/vuthuyduong/fungiclassifiers/tree/master/classification)

##Evaluation

[evaluation](https://github.com/vuthuyduong/fungiclassifiers/tree/master/evaluation)

## Data
There are two datasets available as inputs for fungiclassifiers. The "small" dataset contains ~4000 ITS yeast sequences, checked and validated by the specialists at the Westerdijk Fungal Biodiversity Institute. This dataset were analyzed and released in [Vu D. et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5192050/). The "large" dataset contains ~350K ITS fungal sequences downloaded from GenBank (https://www.ncbi.nlm.nih.gov/) which was used in [Vu D. et al. 2014](https://www.nature.com/articles/srep06837) to evaluate the speed of MLC.

[yeast dataset](http://www.westerdijkinstitute.nl/Download/SmallDatasetOf4KYeastITSSequences.zip) the small demo dataset.

[top50 UNITE](https://github.com/vuthuyduong/fungiclassifiers/tree/master/data/top50UNITE) the large demo dataset. 

## Classifiers

After clustering the DNA sequences by fMLC, the groupings of the sequences can be saved as output of fMLC. A sparse (or complete) similarity matrix (in .sim format) can be saved in the folder where the dataset is given, to capture the similarity structure of the sequences. Based on this similarity matrix, the coordiates of the sequences can be computed and saved (in .outLargeVis format) using LargeVis. Finally, a json file containing the coordinates and metadata of the sequences is resided in the folder DiVE/data folder as an input of DiVE to visualize the data. This json file can be used for visualization by external applications as well.The clustering and visualization results of the two datasets can be found at https://github.com/FastMLC/fMLC/tree/master/data.

## Contact person 

Duong Vu (d.vu@wi.knaw.nl)


