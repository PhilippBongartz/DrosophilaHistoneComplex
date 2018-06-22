#!/usr/bin/env python
# coding: latin1

import random
import time
import sys
import os
import operator
import pickle
import cPickle
import subprocess
import datetime
import numpy as np

######
# 
# This script creates a dataset of reads sampled from a simulated repeat complex.
# The purpose of this dataset is to enable readers of the histone assembly paper
# to run the entire assembly pipeline and assess its general soundness. 
# It is not meant to create a biologically accurate simulation of a repeat complex,
# and for that reason it is probably not useful for quantitative benchmarking. 
# 

def RandomSequence(length):
	seq=''
	for t in range(length):
		seq+='acgt'[int(4*random.random())]
	return seq

def RepeatCopies(seq,copynumber,SNPnumber):
	Copies=[seq for t in range(copynumber)]
	Positions=[int(random.random()*len(seq)) for t in range(SNPnumber)]
	Positions=sorted(Positions)[::-1]
	for t in range(SNPnumber):
		position=Positions[t]
		random.shuffle(Copies)
		rand=int(random.random()*len(Copies)) # Every SNP on another subset
		for x in range(rand):
			Copies[x]=Copies[x][:position]+NotBase[Copies[x][position]][rand%3]+Copies[x][position+1:]
	return Copies		


# Sub 1.4, Ins 11.5, Del 3.4, Match 83.7
NotBase={'a':['c','g','t'],'c':['a','g','t'],'g':['c','a','t'],'t':['c','g','a']}
def PacBioError(seq):
	read=''
	for x in range(len(seq)):
		rand=random.random()
		if rand<0.837+0.115:
			read+=seq[x]
		elif rand<0.837+0.115+0.014:
			read+=NotBase[seq[x]][int(random.random()*3)]
		elif rand<0.837+0.115+0.014+0.034:
			pass
		# Insertion separately otherwise there will always be just one insertion:
		rand=random.random()
		while rand<0.103139:    # geometric formula: 1.0/(1.0-0.103139) - 1.0 = 0.115
			read+='acgt'[int(random.random()*4)]
			rand=random.random()
	return read

# Length distribution taken from the Histone Reads:
Lengthshisto=[0, 323, 427, 411, 355, 353, 358, 321, 293, 321, 281, 275, 241, 239, 226, 
	185, 177, 162, 126, 117, 126, 108, 88, 83, 61, 52, 51, 29, 16, 7, 3, 1, 1, 0, 0, 0, 0, 0, 0, 0]

# Sampling reads according to the length distribution of the histone reads from a given complex sequence
def ReadSampling(coverage, Lengthshisto, genome):
	Lengthsprob=[float(Lengthshisto[t])/float(sum(Lengthshisto)) for t in range(len(Lengthshisto))]
	Lengths=[]
	current_coverage=0
	while current_coverage<coverage:
		rand=random.random()
		length=-1
		prob=0
		while prob<rand:
			length+=1
			prob+=Lengthsprob[length]

		length*=1000
		length+=int(random.random()*1000)
		Lengths.append(length)
		current_coverage=float(sum(Lengths))/float(len(genome))
	Reads=[]
	Starts=[]
	for length in Lengths:
		start=int(random.random()*(len(genome)-length))
		Reads.append(PacBioError(genome[start:start+length]))
		Starts.append(start)
	return Reads,Starts

########
#
# MAIN: Four parameters similar to those observed in the histone complex.
#


# Parameters:
repeatlength=5000
copynumber=100
SNPnumber=200
coverage=90

Repeat=RandomSequence(repeatlength)
Copies=RepeatCopies(Repeat,copynumber,SNPnumber)
genome=''.join(Copies)

# Unique flanking sequences are added to the complex:
FlankingLeft=RandomSequence(10000)
FlankingRight=RandomSequence(10000)
genome=FlankingLeft+genome+FlankingRight


Reads,Starts=ReadSampling(coverage, Lengthshisto, genome)

# The repeat template sequence:
f=open('SimulatedTemplate.fasta','w')
f.write(Repeat+'\n')
f.close()

# The simulated reads sampled with pacbio error from the simulated complex
f=open('SimulatedReads.fasta','w')
for read in Reads:
	f.write('>\n')
	for t in range(0,len(read),100):
		f.write(read[t:t+100]+'\n')
f.close()

# The correct placements of the reads as ground truth
f=open('ReadPlacements','w')
for start in Starts:
	f.write(str(start)+'\n')
f.close()









