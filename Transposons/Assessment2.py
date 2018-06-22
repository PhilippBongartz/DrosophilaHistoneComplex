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
# Takes a Transposon-Signaturefile as argument and calculate the stats.
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
		if sum([sortedcounter[x][1] for x in range(len(sortedcounter))])>10 and sortedcounter[0][1]>2*sortedcounter[1][1]:  # and sortedcounter[0][1]>sum([sortedcounter[x][1] for x in range(len(sortedcounter))])*0.5: 
			konsensus+=sortedcounter[0][0]			
		else:
			konsensus+='x'	# No safe konsensus because coverage too low
			#print sortedcounter
	return konsensus

def KonAcc(Gruppe):	
	konsensus=''
	accuracy=[]
	count=0.0
	first=0.0
	for v in range(len(Gruppe[0])):
		counter={'a':0,'c':0,'g':0,'t':0,'-':0,' ':0}
		for sig in Gruppe:
			counter[sig[v]]+=1
		counter[' ']=0
		sortedcounter=sorted(counter.items(),key=operator.itemgetter(1))
		sortedcounter=sortedcounter[::-1]	
		count+=float(sum([sortedcounter[x][1] for x in range(len(sortedcounter))]))
		first+=float(sortedcounter[0][1])
	if count<0.5:
		return 0,0
	return first,count

def SigDiff(sig1,sig2):
	count=0
	for x in range(len(sig1)):
		if sig1[x]!=sig2[x]:
			count+=1
	return count 		

def Splitter(Gruppe):
	Diffs=[]
	for v in range(len(Gruppe[0])):
		counter={'a':0,'c':0,'g':0,'t':0,'-':0,' ':0}
		for sig in Gruppe:
			counter[sig[v]]+=1
		counter[' ']=0
		sortedcounter=sorted(counter.items(),key=operator.itemgetter(1))
		sortedcounter=sortedcounter[::-1]	
		if sortedcounter[0][1]<1.5*sortedcounter[1][1]:
			Diffs.append(v)
	for sig in Gruppe:
		diffsig=''
		for v in Diffs:
			diffsig+=sig[v]
		print diffsig	
	print 			

# Hier baue ich eine Funktion, die correction accuracy und co berechnet:
# accuracy vs original konsensus, accuracy vs konsensus, F1 vs original konsensus, F1 vs konsensus,
# Und das für small base, big base, all base
# Außerdem collapsed var percentage.
					



SigGroups=[]
def Stats(baxx,counter,Groups,internal):  # Groups kann SigGroups, CorrGroups oder CorrFPGroups sein.
	for t in range(len(counter)):
		if counter[t]>1:
			break
	lower=15
	upper=50		
	# Als erstes stelle ich die big vars und small vars fest:
	countingvars=[{'a':0,'c':0,'g':0,'t':0,'-':0,' ':0} for z in range(len(SigGroups[t][0]))]
	for x in range(len(SigGroups)):
		for y in range(len(SigGroups[x])):
			for z in range(len(SigGroups[x][y])):
				countingvars[z][SigGroups[x][y][z]]+=1
	for z in range(len(SigGroups[t][0])):
		countingvars[z][' ']=0
		#print countingvars[z]

	# Dann mache ich dasselbe für die korrigierten Sigs nur für die Gruppen, die ich auch untersuche: 
	# Warum ist das ok? Weil big var und small var von den Trainingsdaten abhängt
	# Aber die Validierung sollte nur auf den sicheren FlankingClustern passieren.
	relcountingvars=[{'a':0,'c':0,'g':0,'t':0,'-':0,' ':0} for z in range(len(Groups[t][0]))]
	for x in range(len(Groups)):
		for y in range(len(Groups[x])):
			if lower<counter[x]<upper:
				for z in range(len(Groups[x][y])):
					relcountingvars[z][Groups[x][y][z]]+=1
	for z in range(len(Groups[t][0])):
		relcountingvars[z][' ']=0		

	sigaccfirst=0.0
	sigacccount=0.0

	allbigcorrect=0
	allbigcount=0
	allsmallcorrect=0
	allsmallcount=0
	allbasecorrect=0
	allbasecount=0		

	collapsedcount=0
	allcount=0

	smallvarcollection=[[] for z in range(len(SigGroups[t][0]))]

	collapsedbases={'a':0,'c':0,'g':0,'t':0,'-':0,' ':0}
	collapsedto={'a':0,'c':0,'g':0,'t':0,'-':0,' ':0}
	collapsedpos=[0 for jjj in range(len(SigGroups[t][0]))]

	for x in range(len(counter)):
		#if baxx/2<counter[x]<2*baxx:
		if lower<counter[x]<upper:
			bigcorrect=0
			bigcount=0
			smallcorrect=0
			smallcount=0
			basecorrect=0
			basecount=0

			sigkon=Konsensus(SigGroups[x])  # Das ist der Original Konsensus.
			konkon=Konsensus(Groups[x])  

			first,count=KonAcc(SigGroups[x])
			sigaccfirst+=first
			sigacccount+=count

			if internal:
				sigkon=konkon

			for z in range(len(sigkon)):
				allcount+=1
				if sigkon[z]!=konkon[z] and sigkon[z]!='x' and konkon[z]!='x':
					collapsedcount+=1
					collapsedbases[sigkon[z]]+=1
					collapsedto[konkon[z]]+=1
					collapsedpos[z]+=1

					# koncounter={'a':0,'c':0,'g':0,'t':0,'-':0,' ':0,'x':0}
					# for sig in SigGroups[x]:
					# 	koncounter[sig[z]]+=1
					# print koncounter	




			for y in range(len(Groups[x])):
				for z in range(len(Groups[x][y])):
					if Groups[x][y][z]!=' ' and sigkon[z]!='x':   # sigkon=='x' means not enough coverage in Groups[x] for a konsensus.
						basecount+=1
						allbasecount+=1
						
						#if countingvars[z][sigkon[z]]>sum(countingvars[z].values())/3:  # big var
						if countingvars[z][sigkon[z]]==max(countingvars[z].values()):    # majority var
							bigcount+=1
							allbigcount+=1
							if sigkon[z]==Groups[x][y][z]:
								basecorrect+=1
								allbasecorrect+=1
								bigcorrect+=1
								allbigcorrect+=1


						else: 	 # small var  # minority var
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

			# print "Gr{}: Acc: {}, Small: {}, Big: {}.".format(x,float(basecorrect)/float(basecount),
			# 	float(smallcorrect)/float(smallcount),float(bigcorrect)/float(bigcount))


	# Hier berechne ich allsmall davon kann ich dann die allsmallcorrect abziehen und kriege die false positives
	# Für die Berechnung brauche ich eine Sammlung der smallvars die in den Gruppenkonsensen vorkommen.
	allsmall=0
	for z in range(len(SigGroups[t][0])):
		smallvarcollection[z]=list(set(smallvarcollection[z]))
		for base in smallvarcollection[z]:
			allsmall+=relcountingvars[z][base]



	if allbasecount!=0 and allsmallcount!=0 and allbigcount!=0 and allcount!=0:
		print "Acc: {}, Small: {}, Big: {}, collapsed: {}.".format(float(allbasecorrect)/float(allbasecount),
			float(allsmallcorrect)/float(allsmallcount),float(allbigcorrect)/float(allbigcount),float(collapsedcount)/float(allcount))		

	print "bigcount {}, smallcount {}, collapsedcount {}".format(allbigcount,allsmallcount, collapsedcount)
	print "bigcorrect {}, smallcorrect {}".format(allbigcorrect,allsmallcorrect)
	print "Precision small vars: {} = {}/{}".format(float(allsmallcorrect)/float(allsmall),allsmallcorrect,allsmall)
	print "F1-measure for good measure: {}".format(2.0*(float(allsmallcorrect)/float(allsmall))*(float(allsmallcorrect)/float(allsmallcount))/((float(allsmallcorrect)/float(allsmall))+(float(allsmallcorrect)/float(allsmallcount))))
	print "collapsedbases: {}".format(collapsedbases)
	print "collapsedto: {}".format(collapsedto)
	#print "collapsedpos: {}".format(collapsedpos)
	print "SigAccuracy: {}".format(sigaccfirst/sigacccount)





bigsigdiff=0.0
bigcorrdiff=0.0
bigcrossdiff=0.0
bigcount=0.0

# Hier ein kleiner Einschub, der sich die MAs anschaut: LargeScaleVars ? Sind die kleine Cluster unähnlicher ?
# Schaue ich dafür die induzierte MA oder die realignte MA an? Realignt. Konsensus pro Column, dann Alignscore zum Konsensus.
def Alignscore(line,Konsensus):
	count=0
	for x in range(len(line)):
		if Konsensus[x]!=' ' and line[x]!=Konsensus[x]:
			count+=1
	return count


if len(sys.argv)>1 and os.path.exists(sys.argv[1]):
	VarPATH=sys.argv[1]
	liste=VarPATH.split('_')

	try:
		numero=int(liste[-1])
		print "numero:",
		print numero

	except ValueError:
		numero=-1
		print "numero:",
		print numero
		sys.exit()		

	
	if numero<22 and numero>-1:	
		# Dann einlesen und Flankininfo auch einlesen und die Accuracy assessen			
		if (os.path.exists('Signatures_'+str(numero)+'_corr') and os.path.exists('Signatures_'+str(numero)) 
			and os.path.exists('TransposonCopies_'+str(numero)) and os.path.exists('Signatures_'+str(numero)+'_FPcorr')):
			
			f=open('Signatures_'+str(numero),'r')
			Signatures=[]
			for line in f.readlines():
				Signatures.append(line[:len(line)-1])
			f.close()	

			f=open('Signatures_'+str(numero)+'_corr','r')
			SignaturesCorr=[]
			for line in f.readlines():
				SignaturesCorr.append(line[:len(line)-1])
			f.close()	

			f=open('Signatures_'+str(numero)+'_FPcorr','r')
			SignaturesFPCorr=[]
			for line in f.readlines():
				SignaturesFPCorr.append(line[:len(line)-1])
			f.close()

			f=open('TransposonCopies_'+str(numero),'r')
			Flanking=[]
			for line in f.readlines():
				if len(line[:len(line)-1])>0:
					Flanking.append(int(line[:len(line)-1]))
				else:	
					Flanking.append(-1)
			f.close()	

			print len(Signatures),
			print len(SignaturesCorr),
			print len(SignaturesFPCorr),
			print len(Flanking)

			counter=[0 for x in range(max(Flanking)+5)]
			SigGroups=[[] for x in range(max(Flanking)+5)]
			CorrGroups=[[] for x in range(max(Flanking)+5)]
			CorrFPGroups=[[] for x in range(max(Flanking)+5)]
			for t in range(len(Flanking)):
				x=Flanking[t]
				counter[x]+=1
				SigGroups[x].append(Signatures[t])
				CorrGroups[x].append(SignaturesCorr[t])
				CorrFPGroups[x].append(SignaturesFPCorr[t])
			print counter

			# Hier berechne ich mal die normale Coverage:
			zehner=[0 for x in range(20)]
			dreier=[0 for x in range(60)]
			for x in range(len(counter)):
				if counter[x]>-1:
					print min(19,counter[x]/10)
					zehner[min(19,counter[x]/10)]+=1
					dreier[min(59,counter[x]/3)]+=1
			zehner[19]=0
			dreier[59]=0
			maxi=0
			maxx=0
			for x in range(20):
				if zehner[x]>maxi:
					maxi=zehner[x]
					maxx=x
			baxi=0
			baxx=0		
			for x in range(maxx*10,maxx*10+10):
				if dreier[x/3]>baxi:
					baxi=dreier[x/3]
					baxx=x+1
			print zehner		
			print dreier
			print baxx			

			print
			print 
			print "Dataset {}".format(numero)

			print "Stats für Original:"
			Stats(baxx,counter,SigGroups,0)

			print "Stats für CorrFP:"
			Stats(baxx,counter,CorrFPGroups,0)

			print "Stats für CorrFP + internal:"
			Stats(baxx,counter,CorrFPGroups,1)

			print "Stats für Corr:"
			Stats(baxx,counter,CorrGroups,0)

			print "Stats für Corr + internal:"
			Stats(baxx,counter,CorrGroups,1)
			print
			print

			sys.exit()

			


