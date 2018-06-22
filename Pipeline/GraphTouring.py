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
# In this script we do a touring through the graph that is defined by the clusters
# created by Clustering.py. 




Basen2Index={'a':0,'c':1,'g':2,'t':3,'-':4,'_':4,' ':5}
Index2Basen=['a','c','g','t','-',' ']



if len(sys.argv)<2:
  print "Usage: python GraphTouring.py Signatures SeqClass ReadSeqInfo"

#if len(sys.argv)>1 and os.path.exists(sys.argv[1]):

# loading signatures and read information
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

# Loading the ground truth for assessment of the results.
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

  print "Seq:",
  print len(Seqs),
  print "Signatures:",
  print len(CorrectedSignatures)
  print "Reads:",
  print len(Reads)

  # creating ground truth groups
  GroundTruthGruppen=[[] for t in range(102)]
  Cuts=[10000+5000*t for t in range(101)]
  def Place(placement,Cuts,seqlength):
    t=0
    while t<len(Cuts)-1 and placement>Cuts[t]:
      t+=1

    if placement+1000>Cuts[t]:  # determining the cut
      if seqlength<2000:
        pass 
      else:
        t+=1
    return t

  SeqReads=[]
  count=0
  for x in range(len(IndexReads)):
    SeqReads.append([])
    for y in range(len(IndexReads[x])):
      SeqReads[-1].append(count)
      count+=1

  GroundReads=[]
  for x in range(len(IndexReads)):
    GroundReads.append([])
    p=Place(Placements[x],Cuts,len(Seqs[SeqReads[x][0]]))
    for y in range(len(IndexReads[x])):
      GroundTruthGruppen[p+y].append(IndexReads[x][y])
      GroundReads[-1].append(p+y)



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
#   return sum([1 for t in range(len(sig1)) if sig1[t]!=sig2[t] and sig1[t]!=' ' and sig2[t]!=' ']) 

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

ExtReadSeq=[]
for x in range(len(IndexReads)):
  for y in range(len(IndexReads[x])):
    ExtReadSeq.append((x,y))

# We choose initial centroids that don't have lopsided coverage. 
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


# the graph touring:

Clusters=pickle.load(open("Clusters","rb"))
Centroids=pickle.load(open("Centroids","rb"))

if 0:
  # Here we load the unique sequence clusters. 
  f=open('SeqClusters','r')
  UniqueSeqCluster=[]
  for line in f.readlines():
    liste=line.split()
    liste2=[int(l) for l in liste]
    # print liste2
    # for z in liste2:
    #   (x,y)=ExtReadSeq[z]
    #   print KlassReads[x][y],
    # print 
    # for z in liste2:
    #   (x,y)=ExtReadSeq[z]
    #   print y,
    # print 
    # print float(sum([ExtReadSeq[z][1] for z in liste2]))/float(len(liste2))
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

print "Number of Clusters",
print len(Clusters)

# Detecting the left start of the complex:
sortedClusters=sorted([(float(sum([ExtReadSeq[z][1] for z in Clusters[c]]))/float(len(Clusters[c])),c) for c in range(len(Clusters)) if len(Clusters[c])>10])
start=sortedClusters[0][1]

print "Starting cluster {}".format(start)
print float(sum([ExtReadSeq[z][1] for z in Clusters[start]]))/float(len(Clusters[start]))
print [KlassReads[x][y] for (x,y) in [ExtReadSeq[z] for z in Clusters[start]]]

if 0:
  # Mal diese Cluster anschauen:
  Groundtruth=[-1 for z in range(len(FESignatures))]
  for g in range(len(GroundTruthGruppen)):
    for z in GroundTruthGruppen[g]:
      Groundtruth[z]=g

if 0:
  for cluster in Clusters:
    for z in cluster:
      (x,y)=ExtReadSeq[z]
      print Groundtruth[IndexReads[x][y]],
    print 
    print 

  sys.exit()

if 0: # This has been superseded by SeqClustering.c
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

  print len(clusterfirst),
  print len(clusterlast)



# Forward edges
AdjMatrix=[]
for c1 in range(len(Clusters)):
  AdjMatrix.append([])
  for c2 in range(len(Clusters)):
    AdjMatrix[-1].append(len([(x,y) for (x,y) in [ExtReadSeq[z] for z in Clusters[c1] if FESignatures[z][a:-a].count(' ')<covco] if (x,y+1) in [ExtReadSeq[z] for z in Clusters[c2] if FESignatures[z][a:-a].count(' ')<covco]]))

# Backward edges
NAdjMatrix=[]
for c1 in range(len(Clusters)):
  NAdjMatrix.append([])
  for c2 in range(len(Clusters)):
    NAdjMatrix[-1].append(AdjMatrix[c2][c1])

for c1 in range(len(Clusters)):
  for c2 in range(len(Clusters)):
    if NAdjMatrix[c1][c2]!=AdjMatrix[c2][c1]:
      print "Backwards edges do not check out."
      sys.exit()

# Two step forward edges
AAdjMatrix=[]
for c1 in range(len(Clusters)):
  AAdjMatrix.append([])
  for c2 in range(len(Clusters)):
    AAdjMatrix[-1].append(len([(x,y) for (x,y) in [ExtReadSeq[z] for z in Clusters[c1] if FESignatures[z][a:-a].count(' ')<covco] if (x,y+2) in [ExtReadSeq[z] for z in Clusters[c2] if FESignatures[z][a:-a].count(' ')<covco]]))

if 0:
  for c1 in range(len(Clusters)):
    print c1
    for z in Clusters[c1]:
      (x,y)=ExtReadSeq[z]
      print GroundReads[x][y],
    print 
    print 
  sys.exit()

  for c1 in range(len(Clusters)):
    for c2 in range(len(Clusters)):
      if AdjMatrix[c1][c2]:
        print AdjMatrix[c1][c2]
        for z in Clusters[c1]:
          (x,y)=ExtReadSeq[z]
          print GroundReads[x][y],
        print
        for z in Clusters[c2]:
          (x,y)=ExtReadSeq[z]
          print GroundReads[x][y],
        print
        print

  sys.exit()

# A simple layered graph drawing algorithm:
Layers=[[] for l in range(3000)]
Used=[0 for c in range(len(Clusters))]
Layer=[-1 for c in range(len(Clusters))]

Penalty=1
Centroidfrac=1
Avalign=0
Doppel=2

def Scoring(l,c):
  score=0
  p=0
  diff=0
  d=0
  if Layers[l]!=[] or Layers[l-1]!=[] or Layers[min(len(Layers)-1,l+1)]!=[]:
    for cc in range(len(Clusters)):
      if Used[cc]:
        # single step edges
        if AdjMatrix[c][cc]>0:
          if Layer[cc]==l+1:
            score+=AdjMatrix[c][cc]
          else:
            score-=Penalty*AdjMatrix[c][cc]  
            p-=Penalty*AdjMatrix[c][cc]  
        if AdjMatrix[cc][c]>0:
          if Layer[cc]==l-1:
            score+=AdjMatrix[cc][c]
          else:
            score-=Penalty*AdjMatrix[cc][c]
            p-=Penalty*AdjMatrix[cc][c]

        # double step edges
        if AAdjMatrix[c][cc]>0:
          if Layer[cc]==l+2:
            score+=AAdjMatrix[c][cc]*Doppel
            d+=AAdjMatrix[c][cc]*Doppel
          else:
            score-=Penalty*AAdjMatrix[c][cc]*Doppel 
            p-=Penalty*AAdjMatrix[c][cc]*Doppel 
            d-=Penalty*AAdjMatrix[c][cc]*Doppel 
        if AAdjMatrix[cc][c]>0:
          if Layer[cc]==l-2:
            score+=AAdjMatrix[cc][c]*Doppel
            d+=AAdjMatrix[cc][c]*Doppel
          else:
            score-=Penalty*AAdjMatrix[cc][c]*Doppel
            p-=Penalty*AAdjMatrix[cc][c]*Doppel   
            d-=Penalty*AAdjMatrix[cc][c]*Doppel 
            
        if Layer[cc]==l:  # centroid difference
          score-=(Diff(Centroids[cc][a:-a],Centroids[c][a:-a]))/Centroidfrac
          diff+=Diff(Centroids[cc][a:-a],Centroids[c][a:-a])
          p-=(Diff(Centroids[cc][a:-a],Centroids[c][a:-a]))/Centroidfrac

  if len(Layers[l])>0:
    diff/=float(len(Layers[l]))
  if score-p<1.0:  # There have to be edges, just similar centroids is not enough
    return -1.0,-1.0,-1

  return score,p,d  #diff

def MaxGroundTruth(c):
  liste=[GroundReads[ExtReadSeq[z][0]][ExtReadSeq[z][1]] for z in Clusters[c]]
  return (sorted([(liste.count(x),x) for x in liste])[-1])[1]

# The beginning into the first Layer as starting point
Layers[0].append(start)
Used[start]=1
Layer[start]=0

Penalty=1
Centroidfrac=1
Avalign=0


while sum(Used)<len(Clusters):
  maxscore=-10000
  maxdiff=0
  maxc=0
  maxl=0
  maxp=0
  konsc=0
  for c in range(len(Clusters)):
    if not Used[c]:
      for l in range(len(Layers)):
        if 1: #sum([len(Clusters[Layers[l][o]]) for o in range(len(Layers[l]))])+len(Clusters[c])<100:
          score,p,diff=Scoring(l,c)

          if score>maxscore:  # and p==0:
            maxscore=score
            maxl=l
            maxc=c
            maxp=p
            maxdiff=diff



  if maxscore<0:
    print "No more placeable clusters",
    print maxscore
    pickle.dump(Layer,open("Layer","wb"))
    print "Result"
    print Scoring(3,19)
    #sys.exit()
    for l in range(len(Layers)):
      if Layers[l]!=[]:
        print l,
        print len(Layers[l]),
        print Layers[l],
        print [len(Clusters[c]) for c in Layers[l]],
        print [MaxGroundTruth(c) for c in Layers[l]]
    # for l in range(len(Layers)):
    #   if Layers[l]!=[]:
    #     for c in Layers[l]:
    #       for z in Clusters[c]:
    #         (x,y)=ExtReadSeq[z]
    #         print GroundReads[x][y],
    #       print 
    #     print

    sys.exit()

  else:
    print "Score:",
    print maxscore
    print "Cluster {} into Layer {}.".format(maxc,maxl)
    Layers[maxl].append(maxc)
    Used[maxc]=1
    Layer[maxc]=maxl    


print "Result"
print "Layer NumberOfClusters SizeOfClusters GroundTruthGroup"
for l in range(len(Layers)):
  print l,
  print len(Layers[l]),
  print [len(Clusters[c]) for c in Layers[l]],
  print [MaxGroundTruth(c) for c in Layers[l]]
  
pickle.dump(Layer,open("Layer","wb"))
sys.exit()












