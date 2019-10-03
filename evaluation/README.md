# Preprocessing

To evaluate the classification of the models at genus level (classification position = 5) with k-merSize=6, we need to convert the sequences to matrix and split the dataset into training and testing datasets by the following commands:

*fasta2matrix.py 6 yeastITS.fas tmp.matrix*

*addclassification2matrix tmp.matrix yeastITS.classification 5 yeastITS_G.matrix*

*data2trainingandtesting.py yeastITS_G.matrix 4*

Here 4 is the fold number. The training and testing ids will be saved in the results folder.

# Evaluating CNN 

*evaluateCNN.py yeastITS_G.matrix results*/

The evaluation will be saved in the file yeastITS_G.cnn.report

# Evaluating DBN 

*evaluateDBN.py yeastITS_G.matrix results*/

The evaluation will be saved in the file yeastITS_G.dbn.report

# Evaluating BLAST

*evaluateCNN.py yeastITS_G.matrix results*/

The evaluation will be saved in the file yeastITS_G.cnn.report

# Evaluating RDP 

*evaludateCNN.py yeastITS_G.matrix results*/

The evaluation will be saved in the file yeastITS_G.cnn.report


