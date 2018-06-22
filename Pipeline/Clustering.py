#!/usr/bin/env python
# coding: latin1

import random
import numpy as np
import time
import sys
import os
import operator
import cPickle as pickle

#
# In this script we cluster the corrected signatures. 
# 




Basen2Index={'a':0,'c':1,'g':2,'t':3,'-':4,'_':4,' ':5}
Index2Basen=['a','c','g','t','-',' ']



if len(sys.argv)<2:
	print "Usage: python Clustering.py Signatures_corrected SeqClass ReadSeqInfo"

# Loading signatures and read information
if 1:
	SigPath='SimulatedSignatures_corrected'
	SeqClassPath='SimulatedSeqClass'
	ReadSeqPath='SimulatedReadSeqInfo'

	if len(sys.argv)>1 and os.path.exists(sys.argv[1]):
		SigPath=sys.argv[1]
	if len(sys.argv)>2 and os.path.exists(sys.argv[2]):
		SeqClassPath=sys.argv[2]
	if len(sys.argv)>3 and os.path.exists(sys.argv[3]):
		ReadSeqPath=sys.argv[3]  

	CorrectedSignatures=[]
	f=open(SigPath,'r')
	for line in f.readlines():
		if line[-1]=='\n':
			line=line[:len(line)-1]

		CorrectedSignatures.append(line)

	siglength=len(CorrectedSignatures[0])	

	for sig in CorrectedSignatures:
		if len(sig)!=siglength:
			print "Siglengtherror"
			sys.exit()		

	Variations=range(siglength)
	n=siglength-1

	print "siglength: {}".format(siglength)
	f.close()

	f=open(SeqClassPath,'r')
	SeqClass=[]
	for line in f.readlines():
		SeqClass.append(line[:-1])
	f.close()

	f=open(ReadSeqPath,'r')
	IndexReads=[]	
	KlassReads=[]
	for line in f.readlines():
		liste=line.split()
		#print liste

		iiste=[int(l) for l in liste]
		#print iiste

		kiste=[SeqClass[i] for i in iiste]
		#print kiste

		KlassReads.append(kiste)
		IndexReads.append(iiste)

	f.close()		

	count=0
	for x in range(len(IndexReads)):
		new=[]
		for y in range(len(IndexReads[x])):
			if KlassReads[x][y]=='r':
				new.append(count)
				count+=1
			else:
				new.append(-1)
		IndexReads[x]=new

	# We also load the SeqClustering information, to specify the KlassReads: 
	# Here we load the unique sequence clusters. 
	ExtReadSeq=[]
	for x in range(len(IndexReads)):
		for y in range(len(IndexReads[x])):
			ExtReadSeq.append((x,y))

	f=open('SeqClusters','r')
	UniqueSeqCluster=[]
	for line in f.readlines():
	  liste=line.split()
	  liste2=[int(l) for l in liste]
	  UniqueSeqCluster.append(liste2)

	print "Number of UniqueSeqClusters",
	print len(UniqueSeqCluster)

	# The biggest clusters get there own class of unique section:
	sortedNewClusters=sorted([(len(cluster),cluster) for cluster in UniqueSeqCluster])[::-1]
	NewClusters=[cluster for (f,cluster) in sortedNewClusters]
	for c in range(len(NewClusters)):
		if c<7:
			for z in NewClusters[c]:
				(x,y)=ExtReadSeq[z]
				KlassReads[x][y]=['Li','Re','Ze','DA','AN','NN','l'][c]


# Accuracy: 
def Konsensus(liste):
	if liste==[]:
		return ''
	kon=''
	for s in range(len(liste[0])):
		counter=np.array([0,0,0,0,0,0])
		for sig in liste:
			counter[Basen2Index[sig[s]]]+=1
		counter[5]=0
		kon+='acgt- '[np.argmax(counter)]
	return kon

# def Diff(sig1,sig2):
# 	return sum([1 for t in range(len(sig1)) if sig1[t]!=sig2[t] and sig1[t]!=' ' and sig2[t]!=' '])	

mindiff=0
def Diff(sig1,sig2):
	diff=0
	for t in range(len(sig1)):
		if sig1[t]!=' ':
			if sig2[t]!=' ':
				if sig1[t]!=sig2[t]:
					diff+=1
					if diff>mindiff:
						return diff
	return diff

# extending signatures by neighbouring bases
def RealDataSignatureExtension(Signatures, flankingvarnum):
	ExtSignatures=[]
	KlassPattern={'Li':"acgt",'Re':"tgca",'Ze':"aaaa",'DA':"cccc",'AN':"gggg",'NN':"tttt",'l':"ctag"}
	for x in range(len(IndexReads)):
		for y in range(len(IndexReads[x])):
			if KlassReads[x][y] in ['No','Ku','r']:	
				sig=Signatures[IndexReads[x][y]]
			else:
				klasspatt=KlassPattern[KlassReads[x][y]]*siglength
				sig=klasspatt[:siglength]

			linkssig=''
			rechtssig=''
			for z in range(y+1,len(IndexReads[x])):
				if KlassReads[x][z] in ['No','Ku','r']:
					rechtssig+=Signatures[IndexReads[x][z]]
				elif KlassReads[x][z] in KlassPattern:
					klasspatt=KlassPattern[KlassReads[x][z]]*siglength
					rechtssig+=klasspatt[:siglength]					
				else:
					rechtssig+=' '*siglength
			rechtssig+=' '*flankingvarnum		
			for z in range(y-1,-1,-1):
				if KlassReads[x][z] in ['No','Ku','r']:
					linkssig=Signatures[IndexReads[x][z]]+linkssig
				elif KlassReads[x][z] in KlassPattern:
					klasspatt=KlassPattern[KlassReads[x][z]]*siglength
					linkssig=klasspatt[:siglength]+linkssig
				else:
					linkssig=' '*siglength+linkssig
			linkssig=' '*flankingvarnum+linkssig		

			for z in range(flankingvarnum):
				sig=linkssig[len(linkssig)-1-z] + sig + rechtssig[z]	

			ExtSignatures.append(sig)	
	return ExtSignatures

# The corrected signatures are fully extended on both flanks.
FESignatures=RealDataSignatureExtension(CorrectedSignatures, siglength*7)
print "FESignatures have {} sigs, Signatures {}".format(len(FESignatures),len(CorrectedSignatures))
#print set([len(sig) for sig in FESignatures])


# We choose initial centroids that don't have lopsided coverage + large scale vars
weiter=1
a=1
while weiter:
	count=0
	for sig in FESignatures:
		if sig[a]!=' ' and sig[-1*a]!=' ':
			count+=1
	if count*100/len(FESignatures)>10:
		weiter=0
	else:
		a+=1

print "Centroid coverage extent:",
print len(FESignatures[0])-2*a

Centroids=[]
for sig in FESignatures:
  if sig[a]!=' ' and sig[-1*a]!=' ':
    cen=a*' '+sig[a:-a]+a*' '
    if len(cen)!=len(sig):
      print len(cen)
      print len(sig)
      print "!Laengengau"
    Centroids.append(cen)

# Large scale vars don't necessarily have the extented coverage, as is easy to see for the flanking sections:
# So we add a centroid from every large scale var:
z=0
AlreadyUsed=[]
for x in range(len(KlassReads)):
	for y in range(len(KlassReads[x])):
		if KlassReads[x][y]!='r':
			if KlassReads[x][y] not in AlreadyUsed:
				Centroids.append(FESignatures[z])
				AlreadyUsed.append(KlassReads[x][y])
		z+=1

print "Number of centroids",
print len(Centroids)



# The clustering: First round
Clusters=[[] for c in range(len(Centroids))]
for z in range(len(FESignatures)):
  mindiff=10000
  minc=-1
  for c in range(len(Centroids)):
    diff=Diff(FESignatures[z],Centroids[c])
    if diff<mindiff:
      mindiff=diff
      minc=c

  if minc>-1:
    Clusters[minc].append(z)  
    if z%100==0: 	 
      print z,
      sys.stdout.flush()

Centroids2=[]
for cluster in Clusters:
  kon=Konsensus([FESignatures[z] for z in cluster])
  if kon!='':
	  Centroids2.append(kon)


# The clustering: second round
Clusters2=[[] for c in range(len(Centroids2))]
for z in range(len(FESignatures)):
  mindiff=10000
  minc=-1
  for c in range(len(Centroids2)):
    diff=Diff(FESignatures[z],Centroids2[c])
    if diff<mindiff:
      mindiff=diff
      minc=c
  if minc>-1:
    Clusters2[minc].append(z)    
    if z%100==0:
      print z,
      sys.stdout.flush()
print



pickle.dump(Clusters2,open("Clusters","wb"))
pickle.dump(Centroids2,open("Centroids","wb"))

sys.exit()

