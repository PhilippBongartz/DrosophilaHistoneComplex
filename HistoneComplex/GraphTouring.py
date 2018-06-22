#!/usr/bin/env python
# coding: latin1

import string
import operator
import sys
import os
import copy
import time
import pickle
import numpy as np
import itertools
import math
import random
import subprocess

sys.path.append('/Users/Phille/Desktop/Data&Code')
from Functions import *


#
# This script traverses the assembly graph created by the clustering and sorts each cluster into a layer.
# It uses "SavedClustersA2" and creates the file "Layer".
# 

# The histone signatures
VarPATH='UltimateVariations2'
Signatures,siglength,tiefe=SigLoading(VarPATH)

# The corrected signatures
VarPATH="FinalUltimateSigs"
Signatures2,siglength,tiefe=SigLoading(VarPATH)
Variations=range(siglength)

# The information in which read, which signatures and which large-scale vars are contained. 
InfoPATH='NewReadSeqInfointo'
KlassReads,IndexReads=ReadInfo(InfoPATH)


# Seq2ReadSeq=[(-1,-1) for z in range(len(Signatures))]
# for x in range(len(IndexReads)):
#   for y in range(len(IndexReads[x])):
#     if KlassReads[x][y] in ['No','Ku']:
#       Seq2ReadSeq[IndexReads[x][y]]=(x,y)

# if 1:
#   Seq2Read=[-1 for z in range(len(Signatures)+7)]
#   for x in range(len(IndexReads)):
#     for y in range(len(IndexReads[x])):
#       if KlassReads[x][y] in ['Ku','No']:
#         Seq2Read[IndexReads[x][y]]=x


# The MSA
if 0:
  f=open('Histone_MMA')
  MA=[]
  for line in f.readlines():
    MA.append(line[:len(line)-1])
  f.close()

if 1:
  ### InDel-Information is added to the signatures as three vars to have a higher weight.
  weight=3
  for I in range(len(InDels)):
    liste=InDels[I]
    for z in range(len(Signatures2)):
      if z in liste:
        Signatures2[z]+='a'*weight
        Signatures[z]+='a'*weight
      else:
        if 1: #Signatures2[z][InDelCov[I][0]:InDelCov[I][1]].count(' ')==0:
          Signatures2[z]+='t'*weight
          Signatures[z]+='t'*weight
        else:
          Signatures2[z]+=' '*weight
          Signatures[z]+=' '*weight
            
  siglength=len(Signatures2[0]) 


  OldIndexReads=copy.deepcopy(IndexReads)
  OldKlassReads=copy.deepcopy(KlassReads)

  # Unclassified short versions are sorted out.
  count=0
  Kuraus=[]
  for x in range(len(IndexReads)):
    raus=False
    for y in range(len(IndexReads[x])):
      if KlassReads[x][y]=='Ku':
        if IndexReads[x][y] not in atacg and IndexReads[x][y] not in tccta and IndexReads[x][y] not in tcctg:
          raus=True
    if raus:
      IndexReads[x]=[]   
      KlassReads[x]=[]   
      Kuraus.append(x)
      count+=1
  print "{} unclassified short versions are discarded.".format(count)


  ### Here I create some Centroids for the large scale variations. Otherwise they won't feature in the clustering.
  # Centroids are chosen to extend a certain amount of bases into both directions.
  # Large-scale variations don't necessarily have reads with that property,
  # because they might be too long, to accomodate further repeat sequences on both sides.
  # Therefore we combine two reads with a large-scale var to create a centroid:

  # FakeZe:
  IndexReads[4380]+=IndexReads[3684][1:]
  KlassReads[4380]+=KlassReads[3684][1:]


  # 988
  # No No No No AN No
  # +    
  # 948
  # No AN No No No No No

  IndexReads[988]+=IndexReads[948][2:]
  KlassReads[988]+=KlassReads[948][2:]


  # 4865
  # No No No No No DA
  # + 
  # 4557
  # DA No No No NN No No

  IndexReads[4865]+=IndexReads[4557][1:]
  KlassReads[4865]+=KlassReads[4557][1:]

  # 3594
  # Li No No No No Cr
  IndexReads[3594]=[0,0,0,0]+IndexReads[3594]
  KlassReads[3594]=['Li','Li','Li','Li']+KlassReads[3594]

  # 5005
  # No No No No Re
  IndexReads[5005]=IndexReads[5005]+[0,0,0,0]
  KlassReads[5005]=KlassReads[5005]+['Re','Re','Re','Re']

  FESignatures,ExtReadSeq=RealDataSignatureExtension(Signatures2,int(siglength*7),KlassReads,IndexReads)

  print len(FESignatures[0])
  print len(Signatures2[0])

# This parameter specifies the extend of the centroid coverage:
a=1300

count=0
for sig in FESignatures:
  if sig[a]!=' ' and sig[-1*a]!=' ':
    count+=1

print count 

# Loading the clustering results
Clusters=pickle.load(open("SavedClustersA","rb"))
Centroids=pickle.load(open("SavedCentroidsA","rb"))

Clusters2=pickle.load(open("SavedClustersA2","rb"))
Centroids2=[SigKon(cluster,FESignatures) for cluster in Clusters2]    


#############

### reduction of the clusters to long reads with comparable coverage
covco=350  # a coverage cutoff
for c in range(len(Clusters2)):
  Clusters2[c]=[z for z in Clusters2[c] if FESignatures[z][a:-a].count(' ')<covco]
Centroids2=[SigKon(cluster,FESignatures) for cluster in Clusters2]    



##### Touring
# 
# 
#############


#################

# Edge weights of edges between clusters are defined as the number of reads that contain signatures 
# from both clusters in the right order and distance.
# For Forward edges that would be directly adjacent signatures, 
# with the first signature being an element of the first cluster. 
# and the second signature being an element of the second cluster.

# Forward edges
AdjMatrix=[]
for c1 in range(len(Clusters2)):
  AdjMatrix.append([])
  for c2 in range(len(Clusters2)):
    AdjMatrix[-1].append(len([(x,y) for (x,y) in [ExtReadSeq[z] for z in Clusters2[c1] if FESignatures[z][a:-a].count(' ')<covco] if (x,y+1) in [ExtReadSeq[z] for z in Clusters2[c2] if FESignatures[z][a:-a].count(' ')<covco]]))

# Backward edges
NAdjMatrix=[]
for c1 in range(len(Clusters2)):
  NAdjMatrix.append([])
  for c2 in range(len(Clusters2)):
    NAdjMatrix[-1].append(AdjMatrix[c2][c1])

for c1 in range(len(Clusters2)):
  for c2 in range(len(Clusters2)):
    if NAdjMatrix[c1][c2]!=AdjMatrix[c2][c1]:
      print "Falschrum"
      sys.exit()

# Two step forward edges
AAdjMatrix=[]
for c1 in range(len(Clusters2)):
  AAdjMatrix.append([])
  for c2 in range(len(Clusters2)):
    AAdjMatrix[-1].append(len([(x,y) for (x,y) in [ExtReadSeq[z] for z in Clusters2[c1] if FESignatures[z][a:-a].count(' ')<covco] if (x,y+2) in [ExtReadSeq[z] for z in Clusters2[c2] if FESignatures[z][a:-a].count(' ')<covco]]))



# A simple layered graph drawing algorithm:
Layers=[[] for l in range(3000)]
Used=[0 for c in range(len(Clusters2))]
Layer=[-1 for c in range(len(Clusters2))]

Penalty=1   # Penalty for wrong edges
Centroidfrac=1  # Weight of the centroid fit
Avalign=0  # Offset of the centroid fit
Doppel=2  # Score for the longer double step edges.

def Scoring(l,c):
  score=0
  p=0
  diff=0
  d=0
  if Layers[l]!=[] or Layers[l-1]!=[] or Layers[min(len(Layers)-1,l+1)]!=[]:
    for cc in range(len(Clusters2)):
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
          score-=(SigDiff(Centroids2[cc][a:-a],Centroids2[c][a:-a],100)-Avalign)/Centroidfrac
          diff+=SigDiff(Centroids2[cc][a:-a],Centroids2[c][a:-a],100)
          p-=(SigDiff(Centroids2[cc][a:-a],Centroids2[c][a:-a],100)-Avalign)/Centroidfrac

  if len(Layers[l])>0:
    diff/=float(len(Layers[l]))
  if score-p<1.0:  # There have to be edges, just similar centroids is not enough
    return -1.0,-1.0,-1

  return score,p,d  #diff

# For convenience of comparison we start the algorithm at the Layer of the large scale var in the manual assembly.
NN=57
DA=47
NN=51
Ze=38
AN=28
Li=0
Re=112

Start=Li

# We select the starting cluster that corresponds to the large-scale variation.
for c in range(len(Clusters2)):
  if MaxAssAct([ExtReadSeq[z] for z in Clusters2[c] if FESignatures[z][a:-a].count(' ')<covco])==Start:
    Layers[Start].append(c)
    Used[c]=1
    Layer[c]=Start


Penalty=1
Centroidfrac=1
Avalign=0

while sum(Used)<len(Clusters2):
  maxscore=-10000
  maxdiff=0
  maxc=0
  maxl=0
  maxp=0
  konsc=0
  for c in range(len(Clusters2)):
    if not Used[c] and len(Clusters2[c])<90:
      for l in range(len(Layers)):

        if sum([len(Clusters2[Layers[l][o]]) for o in range(len(Layers[l]))])+len(Clusters2[c])<100:
          score,p,diff=Scoring(l,c)

          if score>maxscore:  # and p==0:
            maxscore=score
            maxl=l
            maxc=c
            maxp=p
            maxdiff=diff

  # Calculating the context score for max: This is for diagnostic output only
  vor=[]
  nach=[]
  for cc in range(len(AdjMatrix)):
    if AdjMatrix[maxc][cc]: # and not Used[cc]:
      nach.append((cc,AdjMatrix[maxc][cc]))
    if AdjMatrix[cc][maxc]: # and not Used[cc]:
      vor.append((cc,AdjMatrix[maxc][cc]))

  kontextscore=0.0
  count=0.0
  l=maxl
  if l>1:
    for (cc,adj) in vor:
      sc,pp,dd=Scoring(l-1,cc)
      if sc!=0.0:
        kontextscore+=sc
        count+=1.0
  if l+1<len(Layers):    
    for (cc,adj) in nach:
      sc,pp,dd=Scoring(l+1,cc)
      if sc!=0.0:
        kontextscore+=sc
        count+=1.0

  if count>0.0:  
    kontextscore/=c
  ############################################  

  if maxscore<0:
    print "No more placeable clusters."
    pickle.dump(Layer,open("Layer","wb"))
    sys.exit()

  if 1: #maxscore>0:
    Layers[maxl].append(maxc)
    Used[maxc]=1
    Layer[maxc]=maxl    
    print "Layer {}, Cluster {}, Len {}, Score {}/{}/{}, Ass {}".format(maxl,maxc,len(Clusters2[maxc]),maxscore,maxp,maxdiff,AssAct([ExtReadSeq[z] for z in Clusters2[maxc] if FESignatures[z][a:-a].count(' ')<covco])), 
    print kontextscore
    if len(Clusters2[maxc])>15 and maxl!=MaxAssAct([ExtReadSeq[z] for z in Clusters2[maxc] if FESignatures[z][a:-a].count(' ')<covco]) and MaxAssAct([ExtReadSeq[z] for z in Clusters2[maxc] if FESignatures[z][a:-a].count(' ')<covco])!=-1:
      print "#####################"  
      if 0:  # Diagnostic output:
        # Failure mode:
        for c in range(len(Clusters2)):
          if MaxAssAct([ExtReadSeq[z] for z in Clusters2[c] if FESignatures[z][a:-a].count(' ')<covco])==maxl:
            print c,
            print len(Clusters2[c]),
            print Scoring(maxl,c),
            print MostAssAct([ExtReadSeq[z] for z in Clusters2[c] if FESignatures[z][a:-a].count(' ')<covco],1)
            for cc in range(len(Clusters2)):
              if Used[cc]:
                if AdjMatrix[c][cc] or AdjMatrix[cc][c]:
                  print AdjMatrix[c][cc],
                  print AdjMatrix[cc][c]
                  print Layer[cc],
                  print MostAssAct([ExtReadSeq[z] for z in Clusters2[cc] if FESignatures[z][a:-a].count(' ')<covco],1)
                  print 
            print "Double edge:"      
            for cc in range(len(Clusters2)):
              if Used[cc]:
                if AAdjMatrix[c][cc] or AAdjMatrix[cc][c]:
                  print AAdjMatrix[c][cc],
                  print AAdjMatrix[cc][c]
                  print Layer[cc],
                  print MostAssAct([ExtReadSeq[z] for z in Clusters2[cc] if FESignatures[z][a:-a].count(' ')<covco],1)
                  print

            print
            print       
        sys.exit()


  else:
    for c in range(len(Clusters2)):
      if not Used[c] and sum(NAdjMatrix[c])==0 and len(Clusters2[c])>10:
        maxc=c
        break
    for l in range(len(Layers)):
      if Layers[l]!=[]:
        lasteintrag=l    
    maxl=lasteintrag+50    
    Layers[maxl].append(maxc)
    Used[maxc]=1
    Layer[maxc]=maxl         
    # Output of placement including groundtruth 
    print "Layer {}, Cluster {}, Len {}, Score {}/{}, Groundtruth {}".format(maxl,maxc,len(Clusters2[maxc]),maxscore,maxp,AssAct([ExtReadSeq[z] for z in Clusters2[maxc] if FESignatures[z][a:-a].count(' ')<covco]))  

# And a final output of the resulting assembly including groundtruth. 
for l in range(len(Layers)):
  print "Layer {}:".format(l),
  for c in Layers[l]:
    print "Cluster{},Size{},Groundtruth{}".format(c,len(Clusters2[c]),MaxAssAct([ExtReadSeq[z] for z in Clusters2[c] if FESignatures[z][a:-a].count(' ')<covco])),
  print   


sys.exit()





