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
import datetime


######
# 
# This script calculates the stats for the simulated data set.
#



def Konsensus(Gruppe):	
	konsensus=''
	accuracy=[]
	for v in range(len(Gruppe[0])):
		counter={'a':0,'c':0,'g':0,'t':0,'-':0,' ':0}
		for sig in Gruppe:
			counter[sig[v]]+=1
		counter[' ']=0
		sortedcounter=sorted(counter.items(),key=operator.itemgetter(1))
		sortedcounter=sortedcounter[::-1]	
		konsensus+=sortedcounter[0][0]			
	return konsensus

def SigDiff(sig1,sig2):
	count=0
	for x in range(len(sig1)):
		if sig1[x]!=sig2[x]:
			count+=1
	return count 		

		

# This is a function that calculates correction accuracy and co.
# accuracy vs original konsensus, accuracy vs konsensus, F1 vs original konsensus, F1 vs konsensus,
# And those for small base, big base, all base
# And collapsed var percentage.
# 
if 1:  # Those are simulated data:		
	print "Simulated Data:"			
	f=open('SimulatedCorrVarSigs','r')
	Signatures=[]
	for line in f.readlines():
		Signatures.append(line[:len(line)-1])
	f.close()
	f=open('SimulatedCorrVarSigs_corr','r')
	SignaturesCorr=[]
	for line in f.readlines():
		SignaturesCorr.append(line[:len(line)-1])
	f.close()
	f=open('SimulatedCorrVarSigs_FPcorr','r')
	SignaturesFPCorr=[]
	for line in f.readlines():
		SignaturesFPCorr.append(line[:len(line)-1])
	f.close()

	SigGroups=[[] for x in range(100)]
	CorrGroups=[[] for x in range(100)]
	CorrFPGroups=[[] for x in range(100)]

	for t in range(4999):
		SigGroups[t/50].append(Signatures[t])
		CorrGroups[t/50].append(SignaturesCorr[t])
		CorrFPGroups[t/50].append(SignaturesFPCorr[t])



def Stats(Groups,internal):  # Groups can be SigGroups, CorrGroups or CorrFPGroups.

	counter=[]
	for Group in Groups:
		counter.append(len(Group))

	for t in range(len(counter)):
		if counter[t]>1:
			break	

	# big vars and small vars first
	countingvars=[{'a':0,'c':0,'g':0,'t':0,'-':0,' ':0} for z in range(len(SigGroups[0][0]))]
	for x in range(len(SigGroups)):
		for y in range(len(SigGroups[x])):
			for z in range(len(SigGroups[x][y])):
				countingvars[z][SigGroups[x][y][z]]+=1
	for z in range(len(SigGroups[0][0])):
		countingvars[z][' ']=0
		#print countingvars[z]

	# Now the same for the corrected signatures 
	relcountingvars=[{'a':0,'c':0,'g':0,'t':0,'-':0,' ':0} for z in range(len(Groups[t][0]))]
	for x in range(len(Groups)):
		for y in range(len(Groups[x])):
			if 1:
				for z in range(len(Groups[x][y])):
					relcountingvars[z][Groups[x][y][z]]+=1
	for z in range(len(Groups[t][0])):
		relcountingvars[z][' ']=0		

	allbigcorrect=0
	allbigcount=0
	allsmallcorrect=0
	allsmallcount=0
	allbasecorrect=0
	allbasecount=0		

	collapsedcount=0
	allcount=0

	smallvarcollection=[[] for z in range(len(SigGroups[t][0]))]

	for x in range(len(SigGroups)):
		#if baxx/2<counter[x]<2*baxx:
		if 1:
			bigcorrect=0
			bigcount=0
			smallcorrect=0
			smallcount=0
			basecorrect=0
			basecount=0

			sigkon=Konsensus(SigGroups[x])  # The original consensus.
			konkon=Konsensus(Groups[x])  

			if internal:
				sigkon=konkon

			for z in range(len(sigkon)):
				allcount+=1
				if sigkon[z]!=konkon[z]:
					collapsedcount+=1


			for y in range(len(Groups[x])):
				for z in range(len(Groups[x][y])):
					if Groups[x][y][z]!=' ':
						basecount+=1
						allbasecount+=1
						# if countingvars[z][Groups[x][y][z]]>sum(countingvars[z].values())/3:  # big var
						if countingvars[z][sigkon[z]]>sum(countingvars[z].values())/3:  # big var
							bigcount+=1
							allbigcount+=1
							if sigkon[z]==Groups[x][y][z]:
								basecorrect+=1
								allbasecorrect+=1
								bigcorrect+=1
								allbigcorrect+=1


						else: 	 # small var
							smallvarcollection[z].append(sigkon[z])
							smallcount+=1
							allsmallcount+=1
							if sigkon[z]==Groups[x][y][z]:
								basecorrect+=1
								allbasecorrect+=1
								smallcorrect+=1
								allsmallcorrect+=1

			if bigcount==0:
				bigcount+=1
			if smallcount==0:
				smallcount+=1
			if basecount==0:
				basecount+=1		
			#print "Gr{}: Acc: {}, Small: {}, Big: {}.".format(x,float(basecorrect)/float(basecount),
			#	float(smallcorrect)/float(smallcount),float(bigcorrect)/float(bigcount))

	# Calculating allsmall, substracting allsmallcorrect to get the false positives:
	# First collect the smallvars in the group consensuses
	allsmall=0
	for z in range(len(SigGroups[t][0])):
		smallvarcollection[z]=list(set(smallvarcollection[z]))
		for base in smallvarcollection[z]:
			allsmall+=relcountingvars[z][base]

	print "Acc: {}, Small: {}, Big: {}, collapsed: {}.".format(float(allbasecorrect)/float(allbasecount),
		float(allsmallcorrect)/float(allsmallcount),float(allbigcorrect)/float(allbigcount),float(collapsedcount)/float(allcount))		

	print "bigcount {}, smallcount {}, collapsedcount {}".format(allbigcount,allsmallcount, collapsedcount)
	print "bigcorrect {}, smallcorrect {}".format(allbigcorrect,allsmallcorrect)

	print "Precision small vars: {} = {}/{}".format(float(allsmallcorrect)/float(allsmall),allsmallcorrect,allsmall)
	print "F1-measure for good measure: {}".format(2.0*(float(allsmallcorrect)/float(allsmall))*(float(allsmallcorrect)/float(allsmallcount))/((float(allsmallcorrect)/float(allsmall))+(float(allsmallcorrect)/float(allsmallcount))))


print "Stats Original:"
Stats(SigGroups,0)
print
print
print "Stats CorrFPGroups:"
Stats(CorrFPGroups,0)
print
print
print "Stats CorrGroups:"
Stats(CorrGroups,0)
print
print
print "Stats CorrFPGroups + internal:"
Stats(CorrFPGroups,1)
print
print
print "Stats CorrGroups + internal:"
Stats(CorrGroups,1)
print
print

sys.exit()

