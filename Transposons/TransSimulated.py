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

######
# 
# This script creates a simulated dataset for benchmarking
# Transposons.py --> TransposonReAligner --> TransposonCorrelations.py --> TransNNCorrector.py
# 



def C_IntoAligner(string1, string2):	
	command='./IntoAligner '+string1+' '+string2
	aloutput=subprocess.check_output(command, shell=True)
	info=aloutput.split()
	return info

def MMAker(SeqListe, Vorlage, Pfad, error):
	print "MMAker"
	Columns=[['-' for x in range(len(Vorlage))] for seq in SeqListe] # Columns[seqnumber][vorlagenindex]
	ZwischenColumns=[['' for x in range(len(Vorlage)+1)] for seq in SeqListe]

	seqnumber=-1
	for seq in SeqListe:
		info=C_IntoAligner(Vorlage,seq)  # Wenn der Read länger ist als die Vorlage, dann so rum 
		Editscript=info[1]
		# print Editscript
		Score=info[0]
		if float(Score)/float(len(Vorlage))<error:
			seqnumber+=1
			konindex=0
			seqindex=0
			for x in range(len(Editscript)):
				if Editscript[x]=='m':
					Columns[seqnumber][konindex]=seq[seqindex]		# C_Eintrag[konindex]=string1[seqindex]
					konindex+=1
					seqindex+=1
				if Editscript[x]=='s':	
					Columns[seqnumber][konindex]=seq[seqindex]		# C_Eintrag[konindex]=string1[seqindex]		
					konindex+=1
					seqindex+=1
				
				if Editscript[x]=='d': # i und d werden getauscht, weil Read länger als Vorlage.
					konindex+=1

				if Editscript[x]=='i':		
					ZwischenColumns[seqnumber][konindex]+=seq[seqindex]		# NC_Eintrag[konindex]+=string1[seqindex]
					seqindex+=1
	
	echteseqnumber=seqnumber+1	
	# Jetzt die Ausgabe des MMA: Erstmal die Breite der ZwischenColumns berechnen.
	ZCBreite=[0 for x in range(len(Vorlage)+1)]
	for seqnumber in range(echteseqnumber):
		for y in range(len(Vorlage)+1):
			if len(ZwischenColumns[seqnumber][y])>ZCBreite[y]:
				ZCBreite[y]=len(ZwischenColumns[seqnumber][y])

	string=''
	for seqnumber in range(echteseqnumber):
		string+=ZwischenColumns[seqnumber][0]+(ZCBreite[0]-len(ZwischenColumns[seqnumber][0]))*'-'
		for y in range(len(Vorlage)):
			string+=Columns[seqnumber][y]
			string+=ZwischenColumns[seqnumber][y+1]+(ZCBreite[y+1]-len(ZwischenColumns[seqnumber][y+1]))*'-'
		string+='\n'

	f=open(Pfad,'w')
	f.write(string)
	f.close()	

# Dieses Ding ruft den PW_ReAligner und benennt dann das Ergebnis richtig.
def ReAligner(MMA_string):
	if os.path.exists(MMA_string):
		print "ReAligner:",
		print datetime.datetime.now()
		command='cat '+MMA_string+' | ./PW_ReAligner'
		aloutput=subprocess.check_output(command, shell=True)
		info=aloutput.split()

		f=open('MMA_Real','r')
		g=open('MMA_Simulated_Real','w')
		for line in f.readlines():
			g.write(line)
		f.close()
		g.close()
	else:
		print "Gibt's nicht."


def Correlation(cutoff,MMA_string):
	if not os.path.exists('PyCorrelation'):
		command='gcc -Wall PyCorrelation.c -o PyCorrelation'
		aloutput=subprocess.check_output(command, shell=True)

	if os.path.exists(MMA_string):
		print "Correlation:",
		print "{} for {}mb".format(datetime.datetime.now(),os.stat(MMA_string).st_size/1000000)
		start=time.time()
		command='cat '+MMA_string+' | ./PyCorrelation '+str(cutoff)+' 0 1000000'
		aloutput=subprocess.check_output(command, shell=True)
		info=aloutput.split()
		#print info
		#os.rename('CorrVarSigs','TransSig_'+str(numero))	
		f=open('CorrVarSigs','r')
		g=open('SimulatedCorrVarSigs','w')
		for line in f.readlines():
			g.write(line)
		f.close()
		g.close()	
		print "{} in {} sec".format('SimulatedCorrVarSigs',time.time()-start)

	else:
		print "{} gibt es nicht.".format(MMA_string)	


def Genome_Seq(length):
	f=open("/Users/Phille/Desktop/EcoliDataSets/U00096.2.fasta.txt")
	genome=''
	for line in f.readlines():
		if line[0]!='<':
			genome+=line[:len(line)-1]
	f.close()
	rand=int(random.random()*(len(genome)-length*2))
	return genome[rand:rand+length]

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
		# Insertion separate sonst gibt es immer nur eine Insertion:
		rand=random.random()
		while rand<0.103139:    # geometric formula: 1.0/(1.0-0.103139) - 1.0 = 0.115
			read+='acgt'[int(random.random()*4)]
			rand=random.random()
	return read

# Hier haue ich einfach gelegentlich in 10% 10% mehr Insertions rein:
def VariablePacBioError(seq):
	read=''
	mehr=False
	if random.random()<0.1:
		mehr=True

	for x in range(len(seq)):
		rand=random.random()
		if rand<0.837+0.115:
			read+=seq[x]
		elif rand<0.837+0.115+0.014:
			read+=NotBase[seq[x]][int(random.random()*3)]
		elif rand<0.837+0.115+0.014+0.034:
			pass
		# Insertion separate sonst gibt es immer nur eine Insertion:
		rand=random.random()
		while rand<0.103139:    # geometric formula: 1.0/(1.0-0.103139) - 1.0 = 0.115
			read+='acgt'[int(random.random()*4)]
			rand=random.random()
		if mehr and random.random()<0.10:
			read+='acgt'[int(random.random()*4)]

	return read


def RepeatCopies(seq,copynumber,SNPnumber):
	Copies=[seq for t in range(copynumber)]
	for t in range(SNPnumber):
		position=int(random.random()*len(seq))
		random.shuffle(Copies)
		rand=int(random.random()*len(Copies)) # Every SNP on another subset
		for x in range(rand):
			Copies[x]=Copies[x][:position]+NotBase[Copies[x][position]][rand%3]+Copies[x][position+1:]
	return Copies		

def ReadCreator(Copies,coverage):
	Reads=[]
	for copy in Copies:
		for t in range(coverage):
			Reads.append(PacBioError(copy))
	return Reads		

def VariableReadCreator(Copies,coverage):
	Reads=[]
	for copy in Copies:
		for t in range(coverage):
			Reads.append(VariablePacBioError(copy))
	return Reads

def NNCorrector(VarPATH):
	if os.path.exists(VarPATH):
		print VarPATH
		
		# Dann korrigieren:
		command='python TransNNCorrector.py '+VarPATH
		aloutput=subprocess.check_output(command, shell=True)	



template=Genome_Seq(5000)
print "Template."
Copies=RepeatCopies(template,100,300)
print "Copies."
Reads=ReadCreator(Copies,50)
print "Reads."

Pfad='SimulatedMAcop100cov50'
MMAker(Reads, template, Pfad, 0.50)

ReAligner(Pfad)

MMA_string='MMA_Simulated_Real'

cutoff=10.0
Correlation(cutoff,MMA_string)

VarPfad='SimulatedCorrVarSigs'
NNCorrector(VarPfad)







	


	









