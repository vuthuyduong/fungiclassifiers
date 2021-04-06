#!/usr/bin/env python
# FILE: trainBLAST.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import numpy as np
import os, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='trainBLAST.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile,-p classificationposition -st startingthreshold -et endthreshold -s step",
							   description='''Script that predicts an optimal threshold to separate the sequences based on the given classification''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The output containing a list of threshold and the associately computed F-measure. The optimal threshold is the one producing the highest F-measure.') #optional
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-p','--classificationpos', required=True, type=int, default=0, help='the classification positions for the prediction.')
parser.add_argument('-st','--startingthreshold', type=float, default=0.97, help='starting threshold')
parser.add_argument('-et','--endthreshold', type=float, default=1, help='ending threshold')
parser.add_argument('-s','--step', type=float, default=0.001, help='the step to be increased for the threshold after each step of the prediction.')
parser.add_argument('-mc','--mincoverage', type=int, default=300, help='minimum coverage required for the identitiy of the BLAST comparison.')
args=parser.parse_args()
traindataset= args.input
trainclassificationfilename=args.classification
trainclassificationposition=args.classificationpos
t1 = args.startingthreshold
t2=args.endthreshold
step=args.step
mincoverage=args.mincoverage
modelname=args.out

def GetBase(filename):
	if "." in filename:
		return filename[:-(len(filename)-filename.rindex("."))]
	else:
		return filename

def load_data(modelname,fastafilename,classificationfilename,classificationpos):
	#load seqrecords
	seqids=[]
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load classification
	classificationfile= open(classificationfilename)
	classnames=[]
	newseqrecords=[]
	classes=[]
	for line in classificationfile:
		texts=line.split("\t")
		if line.startswith("#"):
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
			if classname in classnames:
				classes[classnames.index(classname)].append(seqid)
			else:	
				classes.append([seqid])
				classnames.append(classname)
	classificationfile.close()
	basename=GetBase(fastafilename)
	if "/" in basename:
		basename=basename[basename.rindex("/")+1:]
	newfastafilename=modelname + "/" + basename + "." + str(classificationpos) + ".fasta"
	#write selected sequences to file
	SeqIO.write(newseqrecords, newfastafilename, "fasta")
	return newfastafilename,classes



def ComputeSim(fastafilename,seqrecords,mincoverage):
	simmatrix={}
	for seqid in seqrecords.keys():
		simmatrix.setdefault(seqid,{})
		simmatrix[seqid][seqid]=1
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	score=0
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
		i = words[0].rstrip()
		j = words[1].rstrip()
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage	
		if j in simmatrix[i].keys():
			if simmatrix[i][j] < score:
				simmatrix[i][j]=round(score,4)
				simmatrix[j][i]=round(score,4)
		else:
			simmatrix[i][j]=round(score,4)
			simmatrix[j][i]=round(score,4)
		#simmatrix[j][i]=score
	os.system("rm out.txt")
	return simmatrix
		
class PointDef:
	def __init__(self, id,flag,neighbors,seqrecord):
		self.id=id
		self.flag=flag
		self.neighbors = neighbors
		self.seqrecord=seqrecord

class ClusterDef:
	def __init__(self, id, pointids):
		self.id=id
		self.pointids=pointids

def LoadNeighbors(seqids,simmatrix,threshold):
	neighbordict={}
	for seqid in seqids:
		neighbordict.setdefault(seqid, [])
	for i in seqids:
		for j in seqids:
			if not j in simmatrix[i].keys():
				continue
			if simmatrix[i][j] >= threshold:
				if (j not in neighbordict[i]):
					neighbordict[i].append(j)
				if (i not in neighbordict[j]):
					neighbordict[j].append(i)
	#os.system("rm out.txt")
	return neighbordict

def LoadPoints(neigbordict,seqrecords):
	points={}
	for seqid in  seqrecords.keys():
		point = PointDef(seqid,False,neigbordict[seqid],seqrecords[seqid])
		points[seqid]=point
	return points

def ExpandCluster(root,cluster,points):
	for i in root.neighbors:
		if points[i].flag==False:
			points[i].flag=True
			cluster.pointids.append(i)
			ExpandCluster(points[i],cluster,points)		

def Cluster(points,clusters):
	for pointid in points.keys():
		if points[pointid].flag==False:
			points[pointid].flag=True
			cluster = ClusterDef(len(clusters),[])
			cluster.pointids.append(pointid)
			ExpandCluster(points[pointid], cluster,points)
			clusters.append(cluster)
		
def ComputeFmeasure(classes,clusters):
	#compute F-measure
	f=0
	n=0
	for group in classes:
		m = 0
		for cluster in clusters:
			i = len(set(group) & set(cluster.pointids))
			v = float(2*i)/float((len(group) + len(cluster.pointids)))
			if m < v:
				m=v		
		n = n + len(group)
		f = f +	(len(group)*m)	
	return float(f)/float(n) 

def LoadSim(simfilename):
	simmatrix = {} #we use dictionary to reduce the memory constraints 
	simfile = open(simfilename)
	for line in simfile:
		numbers=line.rstrip().split(" ")
		i=numbers[0]
		j=numbers[1]
		if i not in simmatrix.keys():
			simmatrix.setdefault(i, {})
		if j not in simmatrix[i].keys():
			simmatrix[i].setdefault(j, 0)
		if float(numbers[2]) > simmatrix[i][j]:
			simmatrix[i][j]=float(numbers[2])
	simfile.close()		
	return simmatrix

def SaveSim(simmatrix,simfilename):
	simfile=open(simfilename,"w")
	for i in simmatrix.keys():
		for j in simmatrix[i].keys():
			simfile.write(str(i) + " " + str(j) + " " + str(simmatrix[i][j]) + "\n")
	simfile.close()
	
def PredictOpt(predictiondict,simmatrix,records,classes,t1,t2,step,mincoverage):
	optthreshold=0
	bestFmeasure=0
	fmeasuredict={}
	if 'cut-off' in predictiondict.keys():	
		optthreshold=predictiondict['cut-off']
	if 'confidence' in predictiondict.keys():		
		bestFmeasure=predictiondict['confidence']
	if 'fmeasures' in predictiondict.keys():			
		fmeasuredict=predictiondict['fmeasures']
	
	sys.setrecursionlimit(len(records)*2)
	#seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	t=round(t1,4)
	while t < t2:
		print("Compute F-measure for " + str(t) + ":")
		fmeasure=0
		if str(t) in fmeasuredict.keys():
			fmeasure=fmeasuredict[str(t)]
		else:	
			#compute fmeasure
			neighbordict = LoadNeighbors(records.keys(),simmatrix,t)	
			points=LoadPoints(neighbordict,records)
			clusters=[]
			Cluster(points,clusters)	
			fmeasure=round(ComputeFmeasure(classes,clusters),4)
			fmeasuredict[t]=fmeasure
		if fmeasure > bestFmeasure:
			bestFmeasure=fmeasure
			optthreshold=t
		print(str(fmeasure))	
		t=round(t+step,4)
	predictiondict['cut-off']=optthreshold
	predictiondict['confidence']=bestFmeasure
	predictiondict['Number of sequences']=len(records)
	predictiondict['Number of groups']=len(classes)
	predictiondict['fmeasures']=fmeasuredict


def LoadPrediction(optfilename):
	existingprediction={}
	#load classes
	if os.path.exists(optfilename):
		with open(optfilename) as json_file:
			existingprediction = json.load(json_file)
	return existingprediction

def SavePrediction(predictiondict,optfilename):
	#save the whole prediction file
	with open(optfilename,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(predictiondict,json_file)
		else:
			json.dump(predictiondict,json_file,encoding='latin1')
		
def SaveConfig(configfilename,newtraindataset,simfilename,optfilename,classificationfilename,classificationposition,t1,t2,step,mincoverage):
	if not newtraindataset.startswith("/"):
		newtraindataset=os.getcwd() + "/" + newtraindataset
	if not optfilename.startswith("/"):
		optfilename=os.getcwd() + "/" + optfilename
	if not classificationfilename.startswith("/"):
		classificationfilename=os.getcwd() + "/" + classificationfilename
	model="blast"
	#save the config:classifierfilename, model, classificationfilename,classificationpos,k-mer
	configfile=open(configfilename,"w")
	configfile.write("Model: " + model + "\n")
	configfile.write("Optimal threshold prediction filename: " + optfilename + "\n")
	configfile.write("Begin threshold: " + str(t1) + "\n")
	configfile.write("End threshold: " + str(t2) + "\n")
	configfile.write("Step: " + str(step) + "\n")
	configfile.write("Min. coverage: " + str(mincoverage) + "\n")
	configfile.write("Fasta filename: " + newtraindataset + "\n")
	configfile.write("Classification filename: " + classificationfilename + "\n")
	configfile.write("Column number to be classified: " + str(classificationposition) + "\n")
	configfile.write("Similarity filename: " + simfilename + "\n")
	configfile.close()	

##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]
filename=GetBase(traindataset)
if modelname==None or modelname=="":
	modelname=filename.replace(".","_") + "_blast_classifier"
if not os.path.isdir(modelname):
	os.system("mkdir " + modelname)
basename=modelname
if "/" in modelname:
	basename=modelname[modelname.rindex("/")+1:]
optfilename=modelname + "/" + basename + "." + str(trainclassificationposition) +".opt"
simfilename=modelname + "/" + basename + ".sim"

#load data
newtraindataset,classes = load_data(modelname,traindataset,trainclassificationfilename,trainclassificationposition)

#predict optimal threshold	
optthreshold=t1
bestFmeasure=0
if t1 < t2:
	#load sequences
	records=SeqIO.to_dict(SeqIO.parse(newtraindataset, "fasta"))
	#load similarity matrix
	if os.path.exists(simfilename):
		print("Loading similarity matrix " + simfilename)
		simmatrix=LoadSim(simfilename)
	else:	
		print("Computing similarity matrix...")
		simmatrix=ComputeSim(traindataset,records,mincoverage)
		print("Save similarity matrix " + simfilename)
		SaveSim(simmatrix,simfilename)
	predictiondict=LoadPrediction(optfilename)
	PredictOpt(predictiondict,simmatrix,records,classes,t1,t2,step,mincoverage)
	SavePrediction(predictiondict,optfilename)
	#save config	
	configfilename=modelname + "/" + basename + ".config"
	SaveConfig(configfilename,newtraindataset,simfilename,optfilename,trainclassificationfilename,trainclassificationposition,t1,t2,step,mincoverage)
	
	print("The optimal threshold for classification at position " + str(trainclassificationposition) + ": " + str(predictiondict['cut-off']) + ". Best F-measure: " + str(predictiondict['confidence']))
	print("The results is saved in the folder " + modelname + ".") 
	

