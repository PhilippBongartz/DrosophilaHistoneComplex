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



#
# This script detects larger InDels which are added to the corrected signatures before clustering. 
# The detection is based on a preliminary clustering, which greatly speeds up everything,
# but in principle the analysis is possible based on the multiple sequence alignment alone. 
# The idea is to cluster sections of the MSA and analyse the clusters for additional or missing bases. 
#



InfoPATH='NewReadSeqInfointo'
print InfoPATH

f=open(InfoPATH,'r')

KlassReads=[]
IndexReads=[]

for line in f.readlines():
	liste=line.split()
	# print liste
	kiste=[]
	for x in range(3,len(liste)):
		kiste.append(liste[x][:2])
	KlassReads.append(kiste)

	iiste=[]	
	for x in range(3,len(liste)):
		iiste.append(int(liste[x][2:]))
	IndexReads.append(iiste)
f.close()


# Hier die bisherige Assembly:
AlleGruppen=pickle.load(open("pickledAssemblyGruppen","rb"))
AlleGruppen=AlleGruppen[:len(AlleGruppen)-1]

base2index={'a':0,'c':1,'g':2,'t':3,'-':4,' ':5}

f=open('Histone_MMA')
MA=[]
for line in f.readlines():
	MA.append(line[:len(line)-1])
f.close()


def SigDiff(sig1,sig2):
	return len([x for x in range(len(sig1)) if (sig1[x]!=sig2[x] and (sig1[x]!=' ' and sig2[x]!=' '))])

def SegmentKon(Cluster, start, ende):
	kon=''
	for v in range(start,ende):
		counter=[0,0,0,0,0,0]
		for x in Cluster:
			counter[base2index[MA[x][v]]]+=1
		counter[5]=0	
		if max(counter)==0:
			kon+=' '
		else:		
			kon+='acgt-'[counter.index(max(counter))]
	return kon 	

def Clustering(SeqListe,start,ende,K):
	Sigs=[]
	Centroids=[]
	for x in SeqListe:
		sig=''
		for v in range(start,ende):
			sig+=MA[x][v]
		Sigs.append(sig)
		Centroids.append(sig)

	Centroids=list(set(Centroids))	

	# Initializing:
	Clusters=[[] for x in range(len(Centroids))]
	for s in range(len(Sigs)):
		mini=len(sig)
		for x in range(len(Centroids)):
			cen=0
			sd=SigDiff(Centroids[x],Sigs[s])
			if sd<mini:
				mini=sd
				cen=x
		Clusters[cen].append(s)		

	# print [len(cluster) for cluster in Clusters]	
	# print sorted([(len(Clusters[c]),c) for c in range(len(Centroids))])[::-1][:K]

	CClusters=[]
	for (C,s) in sorted([(len(Clusters[c]),c) for c in range(len(Centroids))])[::-1][:K]:
		if C>5:
			CClusters.append(Clusters[s])
	Clusters=CClusters

	Centroids=[]
	for x in range(len(Clusters)):
		Centroids.append(SegmentKon(Clusters[x],start,ende))

	for t in range(1):
		# Clustering:
		Clusters=[[] for x in range(len(Centroids))]
		for s in range(len(Sigs)):
			mini=len(sig)+1
			for x in range(len(Centroids)):
				cen=0
				sd=SigDiff(Centroids[x],Sigs[s])
				if sd<mini:
					mini=sd
					cen=x
			Clusters[cen].append(s)		
		for x in range(len(Centroids)):
			Centroids[x]=SegmentKon(Clusters[x],start,ende)	
	return [[SeqListe[x] for x in Clusters[y]] for y in range(len(Clusters))]		

if 0:  # Here without the preliminary clusters
	for v in range(6500,64200,100):
		w=v+200
		clusters=Clustering(range(len(MA)),v,w,100)
		Konsensen=[SegmentKon(cluster,v,w) for cluster in clusters]
		Sammlungen=[]
		InDelFenster=[]

		for kon in set(Konsensen):
			print kon,
			print kon.count('-')
		print 	


		countsliste=[kon.count('-') for kon in Konsensen if kon.count(' ')==0]
		normalcount=sorted([(Konsensen.count(kon),kon) for kon in Konsensen if kon.count(' ')==0])[-1][0]
		if countsliste!=[] and max(abs(normalcount-max(countsliste)),abs(normalcount-min(countsliste)))>1*((w-v)/200.0):
			print v,
			print countsliste
			konkon=sorted([(Konsensen.count(kon),kon) for kon in Konsensen if kon.count(' ')==0])[-1][1]
			print konkon
			metasammlung=[]
			for kon in Konsensen:
				#if abs(kon.count('-')-konkon.count('-'))>4 and kon.count(' ')==0:
				if SigDiff(konkon, kon)>1*((w-v)/200.0):
					#print kon
					# Hier sammle ich alle Seqs die besser dahin passen und schaue ob die gut sind:
					sammlung=[]
					for z in range(len(MA)):
						segm=MA[z][v:w]
						if SigDiff(segm,konkon)>SigDiff(segm,kon)+1:
							sammlung.append(z)
					#print "Anzahl {}".format(len(sammlung))		
					# for Group in AlleGruppen:
					# 	print len([IndexReads[x][y] for (x,y) in Group if (KlassReads[x][y] in ['Ku','No'] and IndexReads[x][y] in sammlung)]),
					# print 	
					#Sammlungen.append(sammlung)
					metasammlung+=sammlung
			metasammlung=list(set(metasammlung))
			
			if len(metasammlung)>0:
				InDelFenster.append(v)
				# Nochmal ausgeben:
				for kon in Konsensen:
					#if abs(kon.count('-')-konkon.count('-'))>4 and kon.count(' ')==0:
					if SigDiff(konkon, kon)>1*((w-v)/200.0):
						print kon
						# Hier sammle ich alle Seqs die besser dahin passen und schaue ob die gut sind:
						sammlung=[]
						for z in range(len(MA)):
							segm=MA[z][v:w]
							if SigDiff(segm,konkon)>SigDiff(segm,kon)+1:
								sammlung.append(z)
						print "Anzahl {}".format(len(sammlung))		
						for Group in AlleGruppen:
							print len([IndexReads[x][y] for (x,y) in Group if (KlassReads[x][y] in ['Ku','No'] and IndexReads[x][y] in sammlung)]),
						print 			
								
				print "Metaanzahl {}".format(len(metasammlung))		
				Sammlungen.append(metasammlung)
				for Group in AlleGruppen:
					print len([IndexReads[x][y] for (x,y) in Group if (KlassReads[x][y] in ['Ku','No'] and IndexReads[x][y] in metasammlung)]),
				print 		
				print
				print 
			# Hier der Test:
		if not v%1000:
				print v		

		print "Das sind die Sammlungen:"
		print Sammlungen
		print "Das sind die InDelFenster:"
		print InDelFenster



if 1: # This is based on a preliminary clustering of signatures.

	# Hier versuche ich Vars und InDels anhand der cluster zu finden:
	clusters=pickle.load(open('clusters','rb'))
	# Dafür werden sie erstmal in seq übersetzt:
	zclusters=[]
	for cluster in clusters:
		zcluster=[]
		for (x,y) in cluster:
			if KlassReads[x][y] in ['No']: #['Ku','No']:
				zcluster.append(IndexReads[x][y])
		if len(zcluster)>5:
			zclusters.append(zcluster)

	Sammlungen=[]
	InDelFenster=[]
	for v in range(6500,64200,100):
		w=v+200
		# w=v+200
		Konsensen=[]
		for cluster in zclusters:
			if len(cluster)>30:
				reclusters=Clustering(cluster,v,w,len(cluster)/15)
				for recluster in reclusters:
					if len(recluster)>7:
						Konsensen.append(SegmentKon(recluster, v,w))
			else:
				if len(cluster)>7:
					Konsensen.append(SegmentKon(cluster, v,w))	

		countsliste=[kon.count('-') for kon in Konsensen if kon.count(' ')==0]
		normalcount=sorted([(Konsensen.count(kon),kon) for kon in Konsensen if kon.count(' ')==0])[-1][0]
		if countsliste!=[] and max(abs(normalcount-max(countsliste)),abs(normalcount-min(countsliste)))>1*((w-v)/200.0):
			print v,
			print countsliste
			konkon=sorted([(Konsensen.count(kon),kon) for kon in Konsensen if kon.count(' ')==0])[-1][1]
			print konkon
			metasammlung=[]
			for kon in Konsensen:
				#if abs(kon.count('-')-konkon.count('-'))>4 and kon.count(' ')==0:
				if SigDiff(konkon, kon)>1*((w-v)/200.0):
					#print kon
					# Hier sammle ich alle Seqs die besser dahin passen und schaue ob die gut sind:
					sammlung=[]
					for z in range(len(MA)):
						segm=MA[z][v:w]
						if SigDiff(segm,konkon)>SigDiff(segm,kon)+1:
							sammlung.append(z)
					#print "Anzahl {}".format(len(sammlung))		
					# for Group in AlleGruppen:
					# 	print len([IndexReads[x][y] for (x,y) in Group if (KlassReads[x][y] in ['Ku','No'] and IndexReads[x][y] in sammlung)]),
					# print 	
					#Sammlungen.append(sammlung)
					metasammlung+=sammlung
			metasammlung=list(set(metasammlung))
			
			if len(metasammlung)>0:
				InDelFenster.append(v)
				# Nochmal ausgeben:
				for kon in Konsensen:
					#if abs(kon.count('-')-konkon.count('-'))>4 and kon.count(' ')==0:
					if SigDiff(konkon, kon)>1*((w-v)/200.0):
						print kon
						# Hier sammle ich alle Seqs die besser dahin passen und schaue ob die gut sind:
						sammlung=[]
						for z in range(len(MA)):
							segm=MA[z][v:w]
							if SigDiff(segm,konkon)>SigDiff(segm,kon)+1:
								sammlung.append(z)
						print "Anzahl {}".format(len(sammlung))		
						for Group in AlleGruppen:
							print len([IndexReads[x][y] for (x,y) in Group if (KlassReads[x][y] in ['Ku','No'] and IndexReads[x][y] in sammlung)]),
						print 			
								
				print "Metaanzahl {}".format(len(metasammlung))		
				Sammlungen.append(metasammlung)
				for Group in AlleGruppen:
					print len([IndexReads[x][y] for (x,y) in Group if (KlassReads[x][y] in ['Ku','No'] and IndexReads[x][y] in metasammlung)]),
				print 		

				print
				print 

			# Hier der Test:


		if not v%1000:
			print v		

	print "Das sind die Sammlungen:"
	print Sammlungen
	print "Das sind die InDelFenster:"
	print InDelFenster


