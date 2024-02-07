# Fungiclassifiers

# News:
The codes for creating CNN and DBN classifiers have been updated with Tensorflow instead of using Keras.
For classifying the sequences based on BLAST, we have developed a new tool called [dnabarcoder](https://github.com/vuthuyduong/dnabarcoder) that enables us to compute similarity cut-offs for different clades of reference/training datasets. Based on the computed similarity cut-offs we can assign the sequences to the taxonomic group of their best matches.

The code for creating CNN and DBN classifiers has been updated to use TensorFlow instead of Keras. Additionally, for classifying sequences based on BLAST, we have developed a new tool called [dnabarcoder](https://github.com/vuthuyduong/dnabarcoder). This tool enables us to compute similarity cut-offs for different clades of reference/training datasets. Based on these computed similarity cut-offs, we can assign sequences to the taxonomic group of their best matches. It was shown in Vu et al. (2022) that dnabarcoder improves significantly the accuracy and precision of classification.

Fungiclassifiers is the official implementation for the evaluation of the BLAST classification (Altschul et al. 1977), naive Bayesian RDP classifier (Wang et al. 2007), and two deep learning based classifiers using two different models namely convolutional neural network (CNN) (LeCun et al. 2015) and deep belief network (DBN) (Hinton et al. 2006). We aimed to study if deep learning can improve fungal classification. Fungiclassifiers consists of four folders: model, evaluation, data, and Classifiers.

The [models](https://github.com/vuthuyduong/fungiclassifiers/tree/master/models) folder consists of training and classifying scripts of each of the classification models.

The [data](https://github.com/vuthuyduong/fungiclassifiers/tree/master/data) folder consists of the three datasets: the yeast and mould datasets which were checked and validated by the specialists at the Westerdijk Fungal Biodiversity Institute, and the "Top most wanted fungi" (UNITE Community 2017) for the evaluation on different classification methods. The yeast dataset consisting of ~4000 ITS sequences that has been used in (Vu et al. 2018) and released as a subset in (Vu et al. 2016). This dataset was used for the evaluation of different classification methods in which most of the sequences in the test dataset has a label in the train dataset. The mould dataset has been released recently in (Vu et al. 2019). This dataset was to used to train the deep learning models for classifying the "Top most wanted fungi" data whose sequences are unrepresented in the train dataset.

## Citation

Vu, D., Groenewald, M. & Verkley, G. Convolutional neural networks improve fungal classification. Sci Rep 10, 12628 (2020). https://doi.org/10.1038/s41598-020-69245-y

## Pdf

https://www.nature.com/articles/s41598-020-69245-y

## Implementation
The training, classifying in the Classification, and the Evaluation of the BLAST, RDP, CNN and DBN models were implemented in Python 2.7. We used the Keras library (www.keras.io) with tensorflow backend for CNN and the code available at https://github.com/albertbup/deep-belief-network for DBN. The RDP Clas-sifier was downloaded from https://github.com/rdpstaff/classifier.  For the BLAST classification, BLAST version 2.6.0 was installed. 

## Dependencies

Tensorflow

## Contact person 

Duong Vu (d.vu@wi.knaw.nl)

## References
Altschul, S.F., Madden, T.L., Schäffer A.A., et al. (1997). Gapped BLAST and PSI-BLAST: a new generation protein database search programs. Nucleic Acids Res 25, 3389-3402.

Wang, Q., Garrity, G. M., Tiedje, J. M., et al. (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Ap-plied and environmental microbiology 73, 5261–5267. 

LeCun, Y., Bengio, Y., Hinton, G. (2015). Deep learning. Nature, 521, 436-44.

Hinton, G.E., Salakhutdinov, R.R (2006). Reducing the Dimensionality of Data with Neural Networks. Science 313, 504–7.

UNITE Community (2017). UNITE top50 release. Version 01.12.2017. UNITE Community. 

Vu, D., Groenewald, M., Szöke, S., et al. (2016). DNA barcoding analysis of more than 9000 yeast isolates contributes to quantitative thresholds for yeast spe-cies and genera delimitation. Studies in Mycology 85, 91-105. 

Vu, D., Georgievska, S., Szoke, S., et al. (2018). fMLC: Fast Multi-Level Clustering and Visualization of Large Molecular Datasets. Bioinformatics: 1-3.

Vu, D., Groenewald, M., de Vries, M., et al. (2019). Large-scale analysis of fila-mentous fungal DNA barcodes reveals thresholds for species and higher tax-on delimitation. Studies in Mycology 92:135-154.

Vu, D., Nilsson, R.H., Verkley, G. (2022). dnabarcoder: an open-source software package for analyzing and predicting DNA sequence similarity cut-offs for fungal sequence identification. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13651

