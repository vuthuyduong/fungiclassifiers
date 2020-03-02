#!/usr/bin/env python
# FILE: trainingresults2report.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import numpy as np
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

#looking for the best threshold for classification
traindataset=sys.argv[1]
trainclassificationfilename = sys.argv[2]
trainclassificationposition =int(sys.argv[3])
t1=0 #the optimal threshold for classification by clustering in case t2 is not given. Otherwise it is the begin threshold to predict an optimal threshold for clustering
if len(sys.argv) >4: 
	t1 = float(sys.argv[4])
opthreshold=t1
t2=1 #the end threshold used to predict optimal threshold for classification
if len(sys.argv) > 5:
	t2 = float(sys.argv[5])
step = 0.001
if len(sys.argv) > 6:
	step = float(sys.argv[6])
mincoverage=300
if len(sys.argv) > 7:
	mincoverage = int(sys.argv[7])
modelname=""
if len(sys.argv) >8:
	modelname=sys.argv[8]

def load_data(modelname,basename,fastafilename,classificationfilename,classificationposition):
	#load classification
	allseqids=[]
	records= open(classificationfilename)
	allclassification=[]
#	level=""
	for record in records:
		texts=record.split("\t")
#		if record.startswith("#"):
#			if classificationposition < len(texts):
#				level=texts[classificationposition].rstrip()
#			continue 		
		seqid=texts[0].replace(">","").rstrip()
		classname=""
		if classificationposition < len(texts):
			classname=texts[classificationposition].rstrip()
		if classname !="":
			allseqids.append(seqid)
			allclassification.append(classname)
	records.close()
	#load fastafile, save a new fasta file containing only sequences having a classification
	fastafile=open(fastafilename)
	newfastafilename=GetBase(fastafilename)
	if "/" in newfastafilename:
		newfastafilename=newfastafilename[newfastafilename.rindex("/")+1:]
	newfastafilename=modelname + "/" + GetBase(fastafilename) + "."+ str(classificationposition) + ".fasta"
	newfastafile=open(newfastafilename,"w")
	classes=[]
	classnames=[]
	classids=[]
	writetofile=False
	seq=""
	seqid=""
	taxonname=""
	i=0
	for line in fastafile:
		if line.startswith(">"):
			if writetofile==True:
				seqrec=SeqIO.SeqRecord(seq=Seq(seq),id=seqid,description=seqid,name =seqid)
				if taxonname in classnames:
					index=classnames.index(taxonname)
					classes[index].append(seqrec)
					classids[index].append(i)
				else:
					classnames.append(taxonname)
					seqs=[]
					seqs.append(seqrec)
					classes.append(seqs)
					seqids=[]
					seqids.append(i)
					classids.append(seqids)	
			writetofile=False
			seqid=line.split("|")[0].replace(">","").replace(".1","").rstrip()		
			if seqid in allseqids:
				index=allseqids.index(seqid)
				taxonname=allclassification[index]
				if taxonname !="":
					#write to file
					newfastafile.write(">" + seqid + "\n")
					i=i+1
					writetofile=True
		else:
			if writetofile==True:
				seq=line.rstrip()
				newfastafile.write(line)
	if writetofile==True:
		seqrec=SeqIO.SeqRecord(seq=Seq(seq),id=seqid,description=seqid,name =seqid)
		if taxonname in classnames:
			index=classnames.index(taxonname)
			classes[index].append(seqrec)
			classids[index].append(i)
		else:
			classnames.append(taxonname)
			seqs=[]
			seqs.append(seqrec)
			classes.append(seqs)
			seqids=[]
			seqids.append(i)
			classids.append(seqids)			
				
	fastafile.close()
	newfastafile.close()
	return newfastafilename,classes,classnames,classids

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def ComputeBLASTscoreMatrix(fastafilename,seqrecords,mincoverage):
	scorematrix = [[0 for x in range(len(seqrecords))] for y in range(len(seqrecords))] 
	seqids=[]
	for rec in seqrecords:
		seqids.append(rec.id)
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	refid = ""
	score=0
	queryid=""
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
		queryid = words[0].rstrip()
		refid = words[1].rstrip()
		i = seqids.index(queryid)
		j = seqids.index(refid)
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage
		if scorematrix[i][j] < score:
			scorematrix[i][j]=score
			scorematrix[j][i]=score
	os.system("rm out.txt")
	return scorematrix

class PointDef:
	def __init__(self, id,name,flag,neighbors,seqrecord):
		self.id=id
		self.name=name
		self.flag=flag
		self.neighbors = neighbors
		self.seqrecord=seqrecord

class ClusterDef:
	def __init__(self, id, pointids):
		self.id=id
		self.pointids=pointids

def LoadNeighbors(seqrecords,scorematrix,threshold):
	neighborlist=[]
	for i in range(len(seqrecords)):
		neighborlist.append([])
	for i in range(0,len(seqrecords)-1):
		for j in range(i+1,len(seqrecords)):
				if scorematrix[i][j] >= threshold:
					if (j not in neighborlist[i]):
						neighborlist[i].append(j)
					if (i not in neighborlist[j]):
						neighborlist[j].append(i)
	#os.system("rm out.txt")
	return neighborlist

def LoadPoints(neigborlist,seqrecords):
	points=[]
	i=0
	for seqrecord in  seqrecords:
		point = PointDef(i,seqrecord.id,False,neigborlist[i],seqrecord)
		points.append(point)
		i=i+1
	return points

def ExpandCluster(root,cluster,points):
	for i in root.neighbors:
		if points[i].flag==False:
			points[i].flag=True
			cluster.pointids.append(i)
			#print("cluster id:" + str(cluster.id) + "\t" + str(len(cluster.pointids)))
			ExpandCluster(points[i],cluster,points)		

def Cluster(points,clusters):
	for point in points:
		if point.flag==False:
			point.flag=True
			cluster = ClusterDef(len(clusters),[])
			cluster.pointids.append(point.id)
			#print("cluster id:" + str(cluster.id) + "\t" + str(len(cluster.pointids)))
			ExpandCluster(point, cluster,points)
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

def PredictOpt(reportfilename,traindataset,classificationfilename,classificationpos,t1,t2,step,mincoverage,classids):
	#load sequences
	seqrecords=list(SeqIO.parse(traindataset, "fasta"))
	if len(seqrecords)==0:
		return 0,0
	sys.setrecursionlimit(len(seqrecords)*2)
	#seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load similarity matrix
	scorematrix=ComputeBLASTscoreMatrix(traindataset,seqrecords,mincoverage)
	t=t1
	optthreshold=t1
	bestFmeasure=0
		#compute optimal threshold
	optfile=open(reportfilename,"w")
	optfile.write("Threshold\tFmeasure\n")
	while t < t2:
		#load neighbors
		neighborlist = LoadNeighbors(seqrecords,scorematrix,t)	
		points=LoadPoints(neighborlist,seqrecords)
		clusters=[]
		Cluster(points,clusters)	
		fmeasure=ComputeFmeasure(classids,clusters)
		optfile.write(str(t) + "\t" + str(fmeasure) + "\n" )
		if fmeasure > bestFmeasure:
			bestFmeasure=fmeasure
			optthreshold=t
		t=t+step	
	optfile.write("Opt threshold: " + str(optthreshold) + "\tBestFmeasure: " +  str(bestFmeasure) + "\n")
	optfile.close()
	return optthreshold,bestFmeasure

def ComputeVariation(reffilename,mincoverage):
	#load sequeces from the fasta files
	records = list(SeqIO.parse(reffilename, "fasta"))
	scorematrix=ComputeBLASTscoreMatrix(reffilename,records,mincoverage)
	scorelist=[]
	for i in range(0,len(scorematrix)-2):
		for j in range(i+1,len(scorematrix)-1):
			scorelist.append(scorematrix[i][j])
	threshold=1
	minthreshold=1		
	if len(scorelist) >0:
		x=np.array(scorelist)
		minthreshold=min(x)
		threshold=min(x)
		if threshold < np.mean(x)-0.2:
			limit=np.percentile(x,30)
			filtered_x=x[x>=limit]
			threshold=min(filtered_x)
		#threshold=statistics.median(scorelist)
	return threshold,minthreshold

def ComputeVariations(variationfilename,classes,classnames,mincoverage):
	#create json dict
	variations={}
	i=0
	for taxonname in classnames:
		sequences=classes[i]
		if len(sequences) >0:
			fastafilename=taxonname.replace(" ","_") + ".fasta"
			SeqIO.write(sequences,fastafilename,"fasta")
			threshold,minthreshold=ComputeVariation(fastafilename,mincoverage)
			currentvariation={'Threshold': threshold, 'MinThreshold': minthreshold,'NumberOfSequences': len(sequences)}
			variations[taxonname]=currentvariation
			os.system("rm " + fastafilename)
		i=i+1	
	#write to file
	with open(variationfilename,"w") as json_file:
		json.dump(variations,json_file,encoding='latin1')
		
def SaveConfig(configfilename,newtraindataset,optfilename,jsonvariationfilename,classificationfilename,classificationposition,t1,t2,step,mincoverage):
	if not newtraindataset.startswith("/"):
		newtraindataset=os.getcwd() + "/" + newtraindataset
	if not optfilename.startswith("/"):
		optfilename=os.getcwd() + "/" + optfilename
	if not jsonvariationfilename.startswith("/"):
		jsonvariationfilename=os.getcwd() + "/" + jsonvariationfilename
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
	configfile.write("Variation filename: " + jsonvariationfilename + "\n")
	configfile.close()	

##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]
filename=GetBase(traindataset)
if modelname=="":
	modelname=filename.replace(".","_") + "_blast_classifier"
if not os.path.isdir(modelname):
	os.system("mkdir " + modelname)
basename=modelname
if "/" in modelname:
	basename=modelname[modelname.rindex("/")+1:]
optfilename=modelname + "/" + basename + ".opt"

#load data
newtraindataset,classes,classnames,classids = load_data(modelname,basename,traindataset,trainclassificationfilename,trainclassificationposition)

#predict optimal threshold	
optthreshold=t1
if t1 < t2:
	optthreshold,bestFmeasure=PredictOpt(optfilename, newtraindataset,trainclassificationfilename,trainclassificationposition,t1,t2,step,mincoverage,classids)
#compute variation
jsonvariationfilename=modelname + "/" + basename + ".variation"
ComputeVariations(jsonvariationfilename,classes,classnames,mincoverage)
#save config	
configfilename=modelname + "/" + basename + ".config"
SaveConfig(configfilename,newtraindataset,optfilename,jsonvariationfilename,trainclassificationfilename,trainclassificationposition,t1,t2,step,mincoverage)
	
print("The optimal threshold for classification at position " + str(trainclassificationposition) + ": " + str(optthreshold) + ". Best F-measure: " + str(bestFmeasure))
print("The results is saved in the folder " + modelname + ".") 
	

