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
# Here we take a look at the corrected Signatures
# 





Basen2Index={'a':0,'c':1,'g':2,'t':3,'-':4,'_':4,' ':5}
Index2Basen=['a','c','g','t','-',' ']


print
print "Results:"
print 

#if len(sys.argv)>1 and os.path.exists(sys.argv[1]):
if 1:
	SigPath='SimulatedSignatures'
	SeqClassPath='SimulatedSeqClass'
	ReadSeqPath='SimulatedReadSeqInfo'

	Signatures=[]
	f=open(SigPath,'r')
	for line in f.readlines():
		if line[-1]=='\n':
			line=line[:len(line)-1]

		Signatures.append(line)

	siglength=len(Signatures[0])	

	for sig in Signatures:
		if len(sig)!=siglength:
			print "Siglengtherror"
			sys.exit()		

	Variations=range(siglength)
	n=siglength-1

	#print "siglength: {}".format(siglength)
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
	#print "Number of UniqueSeqClusters",
	#print len(UniqueSeqCluster)
	# The biggest clusters get there own class of unique section:
	sortedNewClusters=sorted([(len(cluster),cluster) for cluster in UniqueSeqCluster])[::-1]
	NewClusters=[cluster for (f,cluster) in sortedNewClusters]
	for c in range(len(NewClusters)):
		if c<7:
			for z in NewClusters[c]:
				(x,y)=ExtReadSeq[z]
				KlassReads[x][y]=['Li','Re','Ze','DA','AN','NN','l'][c]

# Here the ground truth information is processed.
if 1:

	PlacementPath='ReadPlacements'
	Placements=[]
	f=open(PlacementPath,'r')
	histo=[0 for t in range(104)]
	for line in f.readlines():
		#print int(line)
		histo[int(line)/5000]+=1
		Placements.append(int(line))
	f.close()

	f=open('SimulatedReads.fasta','r')
	read=''
	Reads=[]
	for line in f.readlines():
		if line[0]=='>' and read!='':
			Reads.append(read)
			read=''
		else:
			read+=line[:-1]
	Reads.append(read)
	f.close()

	f=open('SimulatedSeq.fasta','r')
	read=''
	Seqs=[]
	for line in f.readlines():
		if line[0]=='>' and read!='':
			Seqs.append(read)
			read=''
		else:
			read+=line[:-1]
	Seqs.append(read)
	f.close()

	# Here the signatures are distributed onto ground truth assembly groups:
	# This section still has some magic numbers that would have to be changed,
	# if the Complex was to be simulated with different parameters. 
	GroundTruthGruppen=[[] for t in range(102)]   # 100 copies and 2 flanking sequences
	Cuts=[10000+5000*t for t in range(101)]  # the 101 borders between these sequences
	def Place(placement,Cuts,seqlength):
		t=0
		while t<len(Cuts)-1 and placement>Cuts[t]:
			t+=1

		if placement+1000>Cuts[t]:  # Finding out whether the read was cut 
			if seqlength<2000:
				pass # No
			else:
				t+=1 # Yes
		return t

	SeqReads=[]
	count=0
	for x in range(len(IndexReads)):
		SeqReads.append([])
		for y in range(len(IndexReads[x])):
			SeqReads[-1].append(count)
			count+=1

	# for x in range(len(IndexReads)):
	# 	p=Place(Placements[x],Cuts,len(Seqs[SeqReads[x][0]]))
	# 	for y in range(len(IndexReads[x])):
	# 		GroundTruthGruppen[p+y].append(IndexReads[x][y])

	XY2Truth={}

	GroundReads=[]
	for x in range(len(IndexReads)):
		GroundReads.append([])
		p=Place(Placements[x],Cuts,len(Seqs[SeqReads[x][0]]))
		for y in range(len(IndexReads[x])):
			GroundTruthGruppen[p+y].append(IndexReads[x][y])
			GroundReads[-1].append(p+y)
			XY2Truth[(x,y)]=p+y



CorrectedSignatures=[]
f=open('SimulatedSignatures_corrected','r')
for line in f.readlines():
	CorrectedSignatures.append(line[:-1])
f.close()

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

def Diff(sig1,sig2):
	diff=sum([1 for t in range(len(sig1)) if sig1[t]!=sig2[t] and sig1[t]!=' ' and sig2[t]!=' '])	
	run=sum([1 for t in range(len(sig1)) if sig1[t]!=' ' and sig2[t]!=' '])
	return diff,run

GroundTruthKonsensen=[Konsensus([Signatures[t] for t in GroundTruthGruppen[g] if t>-1]) for g in range(len(GroundTruthGruppen))]
GroundTruthCorKonsensen=[Konsensus([CorrectedSignatures[t] for t in GroundTruthGruppen[g] if t>-1]) for g in range(len(GroundTruthGruppen))]

# Error rate reduction
count=0.0
diffs=0.0
for g in range(1,len(GroundTruthGruppen)-1):
	for t in GroundTruthGruppen[g]:
		if t>-1:
			diff,run=Diff(GroundTruthKonsensen[g],Signatures[t])
			diffs+=diff
			count+=run

print "Error rate of uncorrected Signatures: {} percent.".format((diffs/count)*100)

count=0.0
diffs=0.0
for g in range(1,len(GroundTruthGruppen)-1):
	for t in GroundTruthGruppen[g]:
		if t>-1:
			diff,run=Diff(GroundTruthKonsensen[g],CorrectedSignatures[t])
			count+=siglength
			diffs+=diff
			count+=run

print "Error rate of corrected Signatures: {} percent.".format((diffs/count)*100)
print 

# Internal consistency:
count=0.0
diffs=0.0
for g in range(1,len(GroundTruthGruppen)-1):
	for t in GroundTruthGruppen[g]:
		if t>-1:
			diff,run=Diff(GroundTruthKonsensen[g],Signatures[t])
			diffs+=diff
			count+=run

print "Internal inconsistency of uncorrected Signatures: {} percent.".format((diffs/count)*100)

count=0.0
diffs=0.0
for g in range(1,len(GroundTruthGruppen)-1):
	for t in GroundTruthGruppen[g]:
		if t>-1:
			diff,run=Diff(GroundTruthCorKonsensen[g],CorrectedSignatures[t])
			diffs+=diff
			count+=run

print "Internal inconsistency of corrected Signatures: {} percent.".format((diffs/count)*100)
print 

# small vars error 
#Basen2Index={'a':0,'c':1,'g':2,'t':3,'-':4,' ':5}
Frequency=[[0,0,0,0,0,0] for t in range(siglength)]
for g in range(1,len(GroundTruthGruppen)-1):
	for t in range(len(GroundTruthKonsensen[g])):
		Frequency[t][Basen2Index[GroundTruthKonsensen[g][t]]]+=1

count=0.0
diffs=0.0
for g in range(1,len(GroundTruthGruppen)-1):
	for t in range(siglength):
		if Frequency[t][Basen2Index[GroundTruthKonsensen[g][t]]]<max(Frequency[t]):
			for z in GroundTruthGruppen[g]:
				if z>-1:
					if GroundTruthKonsensen[g][t]!=' ' and ' '!=Signatures[z][t]:
						count+=1.0
						if GroundTruthKonsensen[g][t]!=Signatures[z][t]:
							diffs+=1.0

print "Error rate of small vars in uncorrected Signatures: {} percent.".format((diffs/count)*100)


count=0.0
diffs=0.0
for g in range(1,len(GroundTruthGruppen)-1):
	for t in range(siglength):
		if Frequency[t][Basen2Index[GroundTruthKonsensen[g][t]]]<max(Frequency[t]):
			for z in GroundTruthGruppen[g]:
				if z>-1:
					if GroundTruthKonsensen[g][t]!=' ' and ' '!=CorrectedSignatures[z][t]:
						count+=1.0
						if GroundTruthKonsensen[g][t]!=CorrectedSignatures[z][t]:
							diffs+=1.0

print "Error rate of small vars in corrected Signatures: {} percent.".format((diffs/count)*100)
print 

count=0.0
diffs=0.0
for g in range(1,len(GroundTruthGruppen)-1):
	for t in range(siglength):
		if Frequency[t][Basen2Index[GroundTruthKonsensen[g][t]]]<max(Frequency[t]):
			for z in GroundTruthGruppen[g]:
				if z>-1:
					if GroundTruthKonsensen[g][t]!=' ' and ' '!=Signatures[z][t]:
						count+=1.0
						if GroundTruthKonsensen[g][t]!=Signatures[z][t]:
							diffs+=1.0

print "Internal inconsistency of small vars in uncorrected Signatures: {} percent.".format((diffs/count)*100)


count=0.0
diffs=0.0
for g in range(1,len(GroundTruthGruppen)-1):
	for t in range(siglength):
		if Frequency[t][Basen2Index[GroundTruthKonsensen[g][t]]]<max(Frequency[t]):
			for z in GroundTruthGruppen[g]:
				if z>-1:
					if GroundTruthCorKonsensen[g][t]!=' ' and ' '!=CorrectedSignatures[z][t]:
						count+=1.0
						if GroundTruthCorKonsensen[g][t]!=CorrectedSignatures[z][t]:
							diffs+=1.0

print "Internal inconsistency of small vars in corrected Signatures: {} percent.".format((diffs/count)*100)
print 




# Collapsed variations:
diffs=0
for g in range(1,len(GroundTruthGruppen)-1):
	diffs+=Diff(GroundTruthCorKonsensen[g],GroundTruthKonsensen[g])[0]
print "Of {} group-variation pairs: {} are collapsed.".format((len(GroundTruthGruppen)-2)*siglength,diffs )


# Here we analyse the assembly created by GraphTouring.py

# To do that, we need to transfer the ground truth to the extended signatures
# There are more extended signatures than signatures, because they include large scale variations.
def RealDataSignatureExtension(Signatures, flankingvarnum):
  ExtSignatures=[]
  ExtSig2Sig=[]
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

      if KlassReads[x][y] in ['No','Ku','r']: # sig in extsig as well
        ExtSig2Sig.append(IndexReads[x][y])
      else:  #extsig is not in sig
        ExtSig2Sig.append(-1)

      ExtSignatures.append(sig) 
  return ExtSignatures,ExtSig2Sig

# The corrected signatures are fully extended on both flanks.
FESignatures,ExtSig2Sig=RealDataSignatureExtension(CorrectedSignatures, siglength*7)
#print "FESignatures have {} sigs, Signatures {}".format(len(FESignatures),len(CorrectedSignatures))

Clusters=pickle.load(open("Clusters","rb"))  # The Clusters contain ExtSigs
Centroids=pickle.load(open("Centroids","rb"))  

ExtReadSeq=[]
for x in range(len(IndexReads)):
  for y in range(len(IndexReads[x])):
    ExtReadSeq.append((x,y))

# We determine the area in which centroids are required to have coverage
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

if 0:  # We don't do that anymore, because sequence cluster information is already used in the clustering.
	# Here we load the unique sequence clusters. 
	f=open('SeqClusters','r')
	UniqueSeqCluster=[]
	for line in f.readlines():
	  liste=line.split()
	  liste2=[int(l) for l in liste]
	  print liste2
	  for z in liste2:
	    (x,y)=ExtReadSeq[z]
	    print KlassReads[x][y],
	  print 
	  for z in liste2:
	    (x,y)=ExtReadSeq[z]
	    print y,
	  print 
	  print float(sum([ExtReadSeq[z][1] for z in liste2]))/float(len(liste2))
	  UniqueSeqCluster.append(liste2)

	print len(UniqueSeqCluster)

	# Which of these belong to the left start of the complex? We use those to start the graph touring
	sortedNewClusters=sorted([(float(sum([ExtReadSeq[z][1] for z in cluster]))/float(len(cluster)),cluster) for cluster in UniqueSeqCluster])
	NewClusters=[cluster for (f,cluster) in sortedNewClusters]
	Clusters=NewClusters+Clusters



### reduction of the clusters to long reads with comparable coverage
covco=350
for c in range(len(Clusters)):
	Clusters[c]=[z for z in Clusters[c] if FESignatures[z][a:-a].count(' ')<covco]
Centroids=[Konsensus([FESignatures[z] for z in cluster]) for cluster in Clusters]   

while [] in Clusters:
  Clusters.remove([])
while '' in Centroids:
  Centroids.remove('')

# This has been superseded by the SeqClustering results. 
if 0:
	# The flanking seqs are extracted from the ground truth:
	clusterfirst=[]
	clusterlast=[]
	for z in range(len(ExtReadSeq)):
		(x,y)=ExtReadSeq[z]
		if KlassReads[x][y]=='l':
			if GroundReads[x][y]==0: # y==0 and len(KlassReads[x])>1:
				clusterfirst.append(z)
			if GroundReads[x][y]==101: #y>0 and len(KlassReads[x])-1==y:
				clusterlast.append(z)
	      
	Clusters=[clusterfirst]+Clusters+[clusterlast]
	Centroids=[Konsensus([FESignatures[z] for z in clusterfirst])]+Centroids+[Konsensus([FESignatures[z] for z in clusterlast])]


Layer=pickle.load(open("Layer","rb"))  # Cluster to Layer


# Better ground truth evaluation:
Layers=[[] for t in range(max(Layer)+1)]
for c in range(len(Clusters)):
	for es in Clusters[c]:
		Layers[Layer[c]].append(es)

print 
print "The graph touring assembly, layer by layer:"

correct=0
wrong=0
# We compare each Layers Sigs to the GroundTruth
for l in range(len(Layers)):
	print "Layer {} has".format(l),
	liste=list(set([ExtReadSeq[z] for z in Layers[l]]))
	#print liste[:10]
	liste2=[XY2Truth[(x,y)] for (x,y) in liste]
	#print liste2[:10]

	liste3=sorted([(liste2.count(st),st) for st in set(liste2)])[::-1]
	#print liste2
	print "{} signatures from ground truth group {}, {} from other ground truth groups.".format(liste3[0][0],liste3[0][1],sum([li[0] for li in liste3[1:]]) )
	
	if l==liste3[0][1]:
		correct+=liste3[0][0]
	else:
		wrong+=liste3[0][0]
	wrong+=sum([li[0] for li in liste3[1:]])

	#print liste2
	# for s in Layers[l]:
	# 	print Sig2Truth[s],
	# print 
	# print 

print 
print "{} percent correctly placed signatures.".format( float(int(float(correct*1000)/float(wrong+correct)))/10.0)


sys.exit()





