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

#sys.path.append('/Data&Code')
from Functions import *


#
# This script clusters the extended signatures. 
# To use all available information we first add the detected InDels and the three different versions of 
# the short repeat copy to the signatures. Obviously this is sensible, but it might not be necessary.
# The output of this script is "SavedClustersA2", which is then used in GraphTouring.py
# 

# The Signatures
VarPATH='UltimateVariations2'
Signatures,siglength,tiefe=SigLoading(VarPATH)

# The corrected Signatures
VarPATH="FinalUltimateSigs"
Signatures2,siglength,tiefe=SigLoading(VarPATH)
Variations=range(siglength)

# The information as to which signature and which large-scale variations belongs to which read.
InfoPATH='NewReadSeqInfointo'
KlassReads,IndexReads=ReadInfo(InfoPATH)

# Reads are lists of signatures and large scale vars, this array maps signatures to read and list index.
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
  ### Adding the InDel information
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

# We modify the reads, so here we keep a copy of the original reads.
OldIndexReads=copy.deepcopy(IndexReads)
OldKlassReads=copy.deepcopy(KlassReads)

# Unclassified short versions are sorted out
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
print "{} Ku-Reads aussortiert".format(count)

### Here I add the large scale vars as centroids.
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


# Here we extend signatures by neighbouring bases to improve the resolution of the clustering
# The FullyExtendedSignatures:
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

Centroids=[]
for sig in FESignatures:
  if sig[a]!=' ' and sig[-1*a]!=' ':
    cen=a*' '+sig[a:-a]+a*' '
    if len(cen)!=len(sig):
      print len(cen)
      print len(sig)
      print "!length problem!"
    Centroids.append(cen)

print len(Centroids)

# Now finally the clustering:
# All signatures are distributed among the centroids according to first best fit of the hamming-distance.
if 1:
  Clusters=[[] for c in range(len(Centroids))]
  for z in range(len(FESignatures)):
    mindiff=10000
    minc=-1
    for c in range(len(Centroids)):
      diff=SigDiff(FESignatures[z],Centroids[c],mindiff)
      if diff<mindiff:
        mindiff=diff
        minc=c
    if minc>-1:
      Clusters[minc].append(z)    
      print z,
      sys.stdout.flush()
  print 
  print 
  print 

  # This is an option, where only unique fits are chosen. 
  # Clusters=[[] for c in range(len(Centroids))]
  # for z in range(len(FESignatures)):
  #   mindiff=10000
  #   minc=[-1]
  #   for c in range(len(Centroids)):
  #     diff=SigDiff(FESignatures[z],Centroids[c],mindiff)
  #     if diff==mindiff:
  #       minc.append(c)
  #     elif diff<mindiff:
  #       mindiff=diff
  #       minc=[c]

  #   if len(minc)==1:
  #     Clusters[minc[0]].append(z)    
  #     print z,
  #     sys.stdout.flush()


  #pickle.dump(Clusters,open("SavedClusters","wb"))

  NewClusters=[]
  NewCentroids=[]
  for c in range(len(Clusters)):
    if len(Clusters[c])>3:
      NewClusters.append(Clusters[c])
      NewCentroids.append(Centroids[c])

  # pickle.dump(NewClusters,open("SavedClusters","wb"))
  # pickle.dump(NewCentroids,open("SavedCentroids","wb"))

  pickle.dump(NewClusters,open("SavedClustersA","wb"))
  pickle.dump(NewCentroids,open("SavedCentroidsA","wb"))

# Clusters=pickle.load(open("SavedClusters","rb"))
# Centroids=pickle.load(open("SavedCentroids","rb"))

Clusters=pickle.load(open("SavedClustersA","rb"))
Centroids=pickle.load(open("SavedCentroidsA","rb"))

# Here the consensus of each cluster is calculated as centroids for the second clustering round
Centroids2=[]
for cluster in Clusters:
  kon=SigKon(cluster,FESignatures)
  Centroids2.append(kon)


# This is the second round of clustering:
if 1:
  Clusters2=[[] for c in range(len(Centroids2))]
  for z in range(len(FESignatures)):
    mindiff=10000
    minc=-1
    for c in range(len(Centroids2)):
      diff=SigDiff(FESignatures[z],Centroids2[c],mindiff)
      if diff<mindiff:
        mindiff=diff
        minc=c
    if minc>-1:
      Clusters2[minc].append(z)    
      print z,
      sys.stdout.flush()
  print 

  # The option with unique fits only.
  # Clusters2=[[] for c in range(len(Centroids2))]
  # for z in range(len(FESignatures)):
  #   mindiff=10000
  #   minc=[-1]
  #   for c in range(len(Centroids2)):
  #     diff=SigDiff(FESignatures[z],Centroids2[c],mindiff)
  #     #diff=Passend(FESignatures[z],Centroids2[c])
  #     if diff==mindiff:
  #       minc.append(c)
  #     elif diff<mindiff:
  #       mindiff=diff
  #       minc=[c]
  #   if len(minc)==1:
  #     Clusters2[minc[0]].append(z)    
  #     print z,
  #     sys.stdout.flush()
  # print 


  # The result of the second round of clustering is used in the graph touring. 
  # pickle.dump(Clusters2,open("SavedClusters2","wb"))
  pickle.dump(Clusters2,open("SavedClustersA2","wb"))





