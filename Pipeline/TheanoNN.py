#!/usr/bin/env python
# coding: latin1

import random
import numpy as np
import time
import sys
import os
import operator
import pickle
import cPickle
import theano
import theano.tensor as T
from theano.tensor.nnet import conv
from theano.tensor.nnet import softmax
from theano.tensor import shared_randomstreams


#
# This script is a theano implementation of the correction algorithm. 
# It is the version for the pipeline that processes a simulated repeat complex
# Therefore the KlassReads and IndexReads are adapted for the pipeline outputs.
# 


# Activation functions for neurons
def linear(z): return z
def ReLU(z): return T.maximum(0.0, z)
from theano.tensor.nnet import sigmoid
from theano.tensor import tanh
from theano.tensor.nnet import softmax as rowsoftmax

def softmax(x):
	return rowsoftmax(x.T).T


Basen2Index={'a':0,'c':1,'g':2,'t':3,'-':4,'_':4,' ':5}
Index2Basen=['a','c','g','t','-',' ']

if len(sys.argv)<4:
	print "Usage: python TheanoNN.py Signatures SeqClass ReadSeqInfo threadno"

# Reading in signatures and read information:
if 1:
	SigPath='SimulatedSignatures'
	SeqClassPath='SimulatedSeqClass'
	ReadSeqPath='SimulatedReadSeqInfo'

	if len(sys.argv)>1 and os.path.exists(sys.argv[1]):
		SigPath=sys.argv[1]
	if len(sys.argv)>2 and os.path.exists(sys.argv[2]):
		SeqClassPath=sys.argv[2]
	if len(sys.argv)>3 and os.path.exists(sys.argv[3]):
		ReadSeqPath=sys.argv[3]

	print "SigPath: {}".format(SigPath)
	print "SeqClassPath: {}".format(SeqClassPath)
	print "ReadSeqPath: {}".format(ReadSeqPath)


	Signatures=[]
	f=open(SigPath,'r')
	for line in f.readlines():
		if line[-1]=='\n':
			line=line[:len(line)-1]
		Signatures.append(line)

	siglength=len(Signatures[0])	

	print "signature length: {}".format(siglength)
	f.close()

	count=0
	for sig in Signatures:
		if len(sig)!=siglength:
			print "Siglengtherror",
			print count,
			print len(sig)
			sys.exit()	
		count+=1	

	Variations=range(siglength)
	n=siglength-1


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


# Checking with the ground truth.
if 0:

	PlacementPath='ReadPlacements'
	Placements=[]
	f=open(PlacementPath,'r')
	histo=[0 for t in range(104)]
	for line in f.readlines():
		print int(line)
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
	print len(Signatures)


	print len(Reads)

	print histo
	print min(Placements),
	print max(Placements)

	print len(Placements),
	print len(IndexReads)


	# Sequences are placed into groups according to the placement information
	GroundTruthGruppen=[[] for t in range(102)]
	Cuts=[10000+5000*t for t in range(101)]
	def Place(placement,Cuts,seqlength):
		t=0
		while t<len(Cuts)-1 and placement>Cuts[t]:
			t+=1

		if placement+1000>Cuts[t]: 
			if seqlength<2000:
				pass # Stimmt so
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

	for x in range(len(IndexReads)):
		if Placements[x]<20000:
			print x
			print len(Seqs[SeqReads[x][0]]),
			print len(Seqs[SeqReads[x][-1]])
			print Placements[x]
			print Place(Placements[x],Cuts,len(Seqs[SeqReads[x][0]]))
			print 

	sys.exit()

# A function to extend signatures with neighbouring bases
def RealDataSignatureExtension(Signatures, flankingvarnum):
	ExtSignatures=[]
	KlassPattern={'Li':"acgt",'Re':"tgca",'Ze':"aaaa",'DA':"cccc",'AN':"gggg",'NN':"tttt",'l':"ctag"}
	for x in range(len(IndexReads)):
		for y in range(len(IndexReads[x])):
			if KlassReads[x][y] in ['No','Ku','r']:	

				sig=Signatures[IndexReads[x][y]]
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

#################### Important parameter: the length of the extension
flankingvarnum=siglength
# flankingvarnum=0
length=siglength+flankingvarnum+flankingvarnum
ExtSignatures=RealDataSignatureExtension(Signatures,flankingvarnum)

#Truedata=TrueData(Signatures)
for x in range(len(ExtSignatures)-1):
	if len(ExtSignatures[x])!=len(ExtSignatures[x+1]):
		print len(ExtSignatures[x+1])
print "extended signature length:",
print length





#########################################################
			
# Turns Signatures into vectors
def Data2Vectors(ISignatures):	
	Sigvectoren=[]
	for i in range(len(ISignatures)):
		OneInput=np.zeros((len(ISignatures[0])*5, 1))   # five entries per variation except the target variation
		
		for v in range(len(ISignatures[0])):
			if ISignatures[i][v]!=' ':
				OneInput[Basen2Index[ISignatures[i][v]]+v*5]=1.0
		Sigvectoren.append(OneInput)	
	return Sigvectoren		


# The implementation of the neural network
class Network(object):

	def __init__(self,layersizes, activation_fn1=sigmoid, activation_fn2=sigmoid):
		if len(layersizes)!=3:
			print "Dies ist ein DreierNN."
			sys.exit()

		self.sizes=layersizes
		self.activation_fn1=activation_fn1
		self.activation_fn2=activation_fn2
		self.WeightInitializer()

		# Input, Target and Output are vectors
		self.x = T.dmatrix("x")  
		self.y = T.dmatrix("y")
		
		self.a1 = T.dmatrix("a1") 
		self.a2 = T.dmatrix("a2") 

		self.z0 = T.dmatrix("z0") 
		self.z1 = T.dmatrix("z1")   # This is the output


	def WeightInitializer(self):	

		self.W0 = theano.shared(
			np.asarray(
				np.random.normal(loc=0.0, scale=np.sqrt(1.0/self.sizes[1]), size=(self.sizes[1],self.sizes[0])),
					dtype=theano.config.floatX),
			name='W0', borrow=True)

		self.W1 = theano.shared(
			np.asarray(
				np.random.normal(loc=0.0, scale=np.sqrt(1.0/self.sizes[2]), size=(self.sizes[2],self.sizes[1])),
					dtype=theano.config.floatX),
			name='W1', borrow=True)		

		self.B0 = theano.shared(
			np.asarray(np.random.normal(loc=0.0, scale=1.0,size=(self.sizes[1],1),),
						dtype=theano.config.floatX),
			name='B0', borrow=True)

		self.B1 = theano.shared(
			np.asarray(np.random.normal(loc=0.0, scale=1.0,size=(self.sizes[2],1),),
						dtype=theano.config.floatX),
			name='B1', borrow=True)

		self.params = [self.W0,self.W1,self.B0,self.B1]

	def SGD(self, minibatchsize, trainingsdata, epochs, eta, lambd):
		# The computation graph of the neural network

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		# the regularizazion term
		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd

		act_reg = ((self.z1-0.5)**2).sum()

		# the cost function
		qcost = (((self.z2 - self.y)**2).sum()) * 0.5 + regularization #+ act_reg

		y_out = T.argmax(self.z2, axis=0)
		grads = T.grad(qcost, self.params)

		cost = theano.function([self.x,self.y],qcost)
		feedforward = theano.function([self.x],y_out)
		gradients = theano.function([self.x,self.y],grads)

		changesW0 = T.dmatrix("changesW0")
		updaterW0 = theano.function([changesW0], self.W0, updates=[(self.W0, self.W0+changesW0)])

		changesW1 = T.dmatrix("changesW1")
		updaterW1 = theano.function([changesW1], self.W1, updates=[(self.W1, self.W1+changesW1)])

		changesB0 = T.dmatrix("changesB0")
		updaterB0 = theano.function([changesB0], self.B0, updates=[(self.B0, self.B0+changesB0)])

		changesB1 = T.dmatrix("changesB1")
		updaterB1 = theano.function([changesB1], self.B1, updates=[(self.B1, self.B1+changesB1)])

		# the training
		for i in range(epochs):
			z = random.randint(0,len(trainingsdata)-minibatchsize)

			mini_batch = trainingsdata[z:z+minibatchsize]

			changes = [0.0 for param in self.params] # updating after each minibatch
			for (x,y) in mini_batch:
				upgrads = gradients(x,y)
				changes = [ change - upgrad*eta/float(minibatchsize) for upgrad,change in zip(upgrads,changes)]   # /minibatchsize damit die Lernrate davon unabhängig bleibt.

			#self.params = [ param-change for param, change in zip(self.params, changes) ]
			updaterW0(changes[0])
			updaterW1(changes[1])
			updaterB0(changes[2])
			updaterB1(changes[3])

	# Eperimental code
	def SelftuningSGD(self, minibatchsize, trainingsdata, epochs, eta, lambd, selftuningstep):		
		# Hier wird die Costfunktion in Abhängigkeit von self.x, self.y und den Weights und Biases berechnet.
		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd

		qcost = -T.mean(self.y*T.log(self.z2)+(1-self.y)*T.log(1-self.z2)) + regularization

		starget = T.argmax(self.y,axis=0)
		y_out = T.argmax(self.z2, axis=0)
		grads = T.grad(qcost, self.params)

		cost = theano.function([self.x,self.y],qcost)
		feedforward = theano.function([self.x],y_out)
		gradients = theano.function([self.x,self.y],grads)

		target = theano.function([self.y],starget)

		changesW0 = T.dmatrix("changesW0")
		updaterW0 = theano.function([changesW0], self.W0, updates=[(self.W0, self.W0+changesW0)])

		changesW1 = T.dmatrix("changesW1")
		updaterW1 = theano.function([changesW1], self.W1, updates=[(self.W1, self.W1+changesW1)])

		changesB0 = T.dmatrix("changesB0")
		updaterB0 = theano.function([changesB0], self.B0, updates=[(self.B0, self.B0+changesB0)])

		changesB1 = T.dmatrix("changesB1")
		updaterB1 = theano.function([changesB1], self.B1, updates=[(self.B1, self.B1+changesB1)])

		#hiddenstate = theano.function([self.x],self.z1)

		trainingsdatacorrect=[0 for yyy in range(len(trainingsdata))]	

		for i in range(epochs):
			z = random.randint(0,len(trainingsdata)-minibatchsize)
			mini_batch = trainingsdata[z:z+minibatchsize]

			changes = [0.0 for param in self.params] # Upgedated wird am Ende des Minibatch
			for (x,y) in mini_batch:
				upgrads = gradients(x,y)
				changes = [ change - upgrad*eta/float(minibatchsize) for upgrad,change in zip(upgrads,changes)]  
				# /minibatchsize damit die Lernrate davon unabhängig bleibt.

			#self.params = [ param-change for param, change in zip(self.params, changes) ]
			updaterW0(changes[0])
			updaterW1(changes[1])
			updaterB0(changes[2])
			updaterB1(changes[3])

			# Hier wird eta an die true und false positives angepasst:
			if not i%selftuningstep:
				#avcost=0.0
				correct=0.0
				neuverbessert=0
				neuverschlechtert=0
				for yyy in range(len(trainingsdata)):
					(x,y) = trainingsdata[yyy]
					# print hiddenstate
					# print hiddenstate(x)
					# print

					#avcost+=cost(x,y)

					if feedforward(x)==target(y):
						correct+=1.0
						if not trainingsdatacorrect[yyy]:   # war vorher nicht korrekt korrigiert
							neuverbessert+=1
						trainingsdatacorrect[yyy]=1
					else:		# wurde nicht korrekt korrigiert
						if trainingsdatacorrect[yyy]:  # war aber vorher korrekt korrigiert:
							neuverschlechtert+=1
						trainingsdatacorrect[yyy]=0	

				if neuverschlechtert==0:
					eta*=1.1
				elif 2*neuverschlechtert>neuverbessert:
					eta/=1.1	

				#print "Average Cost = {}.".format(avcost/float(len(trainingsdata)))
				print "Correct: {} prozent. {} neuverbessert, {} neuverschlechtert. Neu-Eta: {}".format(float(int(1000.0*correct/float(len(trainingsdata))))/10.0, neuverbessert, neuverschlechtert,eta)					
				self.Accuracy(trainingsdata,lambd)



	# Outputting accuracy and error
	def Accuracy(self,trainingsdata, lambd):		

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd

		qcost = (((self.z2 - self.y)**2).sum()) * 0.5 + regularization #+ act_reg

		starget = T.argmax(self.y,axis=0)
		y_out = T.argmax(self.z2, axis=0)
		grads = T.grad(qcost, self.params)

		cost = theano.function([self.x,self.y],qcost)
		feedforward = theano.function([self.x],y_out)
		gradients = theano.function([self.x,self.y],grads)
		target = theano.function([self.y],starget)

		output = theano.function([self.x],self.z2)

		avcost=0.0
		correct=0.0

		basehisto=[0]*self.sizes[2]
		truehisto=[0]*self.sizes[2]
		truehits=[0]*self.sizes[2]

		for (x,y) in trainingsdata:
			avcost+=cost(x,y)


			ff=feedforward(x)
			tt=target(y)

			basehisto[int(ff)]+=1
			truehisto[int(tt)]+=1

			if ff==tt:
				truehits[int(tt)]+=1
				correct+=1.0

		print "Average Cost = {}.".format(avcost/float(len(trainingsdata)))
		print "Correct: {} prozent.".format(float(int(1000.0*correct/float(len(trainingsdata))))/10.0)	
		return float(int(1000.0*correct/float(len(trainingsdata))))/10.0

	# predicting a base
	def PredictedBase(self, x):		

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		y_out = T.argmax(self.z2, axis=0)
		
		feedforward = theano.function([self.x],y_out)

		y=feedforward(x)
		Index2Basen=['a','c','g','t','-',' ']

		return Index2Basen[y]

	# predicting many bases
	def PredictedBases(self, InputData):		

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		y_out = T.argmax(self.z2, axis=0)
		
		feedforward = theano.function([self.x],y_out)
		Index2Basen=['a','c','g','t','-',' ']

		NewBases=[]
		Index2Basen=['a','c','g','t','-',' ']

		for x in InputData:
			y=feedforward(x)
			NewBases.append(Index2Basen[int(y)])

		return NewBases
		
	# a diagnostic function
	def AvAct(self,trainingsdata, lambd):		

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd

		qcost = (((self.z2 - self.y)**2).sum()) * 0.5 + regularization #+ act_reg

		starget = T.argmax(self.y,axis=0)
		y_out = T.argmax(self.z2, axis=0)
		grads = T.grad(qcost, self.params)

		cost = theano.function([self.x,self.y],qcost)
		feedforward = theano.function([self.x],y_out)
		gradients = theano.function([self.x,self.y],grads)
		target = theano.function([self.y],starget)

		output = theano.function([self.x],self.z2)

		act1 = theano.function([self.x],self.z1)

		avcost=0.0
		correct=0.0

		basehisto=[0,0,0,0,0]

		avact1=OneInput=np.zeros((self.sizes[1], 1))
		avact2=OneInput=np.zeros((self.sizes[2], 1))

		for (x,y) in trainingsdata:
			avact1 += act1(x)
			avact2 += output(x)

		print "activation1:"
		print avact1/float(len(trainingsdata))
		print "activation2:"
		print avact2/float(len(trainingsdata))



Vektoren=Data2Vectors(ExtSignatures)

# Training of one variation:
def Training(v):
	print
	print
	print "Variation {}:".format(v)
	Var=v+flankingvarnum
	trainingsdata=[(np.concatenate((Vektor[:Var*5],Vektor[(Var+1)*5:]),axis=0),Vektor[Var*5:(Var+1)*5]) for Vektor in Vektoren if np.linalg.norm(Vektor[Var*5:(Var+1)*5])>0.5]   # Var
	namestring=SigPath+str(v)+"_net"
	if os.path.exists("Nets/"+namestring):
		net=cPickle.load(open("Nets/"+namestring,"rb"))
	else:
		net = Network([(length-1)*5,5,5], tanh, sigmoid)
	lambd=1.0 / float((length-1)*5)
	acc=net.Accuracy(trainingsdata, lambd)
	if acc < 97:
		print "Further Training:"
		eta=0.0030000
		net.SGD(10,trainingsdata,500000,eta, lambd)
		net.Accuracy(trainingsdata, lambd)		
	cPickle.dump(net,open("Nets/"+namestring,"w"))

#### Here I try out multithreading, but it doesn't seem to work well with theano
import threading

class myThread (threading.Thread):
	def __init__(self, v):
		threading.Thread.__init__(self)
		self.v = v

	def run(self):
		print "Starting " + self.name
		Training(self.v)
		print "Exiting " + self.name
# try:
# 	thread1 = myThread(0)
# 	thread2 = myThread(1)
# except:
# 	print "Oops"


# print threading.active_count()
# thread1.start()
# print threading.active_count()
# thread2.start()
# print threading.active_count()

################ Here all nets are trained and pickled
if not os.path.exists("Nets"):
	os.makedirs("Nets")

# The threading doesn't really work, we parallelize more primitively.
if 0:
	threadno=1
	if len(sys.argv)>4:
		threadno=int(sys.argv[4])

	if threadno>1:
		v=0
		while v<siglength:
			if threading.active_count()<threadno:
				try:
					thread = myThread(v)
					thread.start()
					v+=1
					print v
				except:
					print "Oops"
				print threading.enumerate()

	else:
		for v in range(0,siglength):  #(125,siglength): #,siglength):
			Training(v)

# More primitive parallelization:
if 1:
	startv=0
	if len(sys.argv)>4:
		startv=int(sys.argv[4])

	# We just go through all variations in turn, starting with startv
	for v in range(startv,siglength+startv):  #(125,siglength): #,siglength):
		v=v%siglength
		Training(v)

################ Here all nets are trained and pickled


################# Corrections:
# Nets are loaded and used to correct the signatures

CorrectedSigListen=[[] for z in range(len(ExtSignatures))] # There are several rounds of correction.



print "Loading nets."
Nets=[]
for v in range(siglength):  #(125,siglength): #,siglength):
	namestring=SigPath+str(v)+"_net"
	if not os.path.exists("Nets/"+namestring):
		print "Nets/"+namestring
		print "does not exist."
		sys.exit()
		
	net=cPickle.load(open("Nets/"+namestring,"rb"))
	Nets.append(net)

print "Correction."
for d in range(5):  # Correction
	print "Round {}.".format(d)
	for z in range(len(ExtSignatures)):
		CorrectedSigListen[z].append('')

	for v in range(siglength):
		#print v,
		sys.stdout.flush()

		Var=v+flankingvarnum
		InputData=[ np.concatenate((Vektor[:Var*5],Vektor[(Var+1)*5:]),axis=0) for Vektor in Vektoren]
		
		NewBases=Nets[v].PredictedBases(InputData)

		for z in range(len(ExtSignatures)):
		
			if Signatures[z][v]==' ':
				CorrectedSigListen[z][-1]+=' '
			else:
				#print z
				sys.stdout.flush()

				CorrectedSigListen[z][-1]+=NewBases[z]
			
	#print		
	# The new signatures are extended and then vectorized:
	ExtSignatures=RealDataSignatureExtension([SigListe[d] for SigListe in CorrectedSigListen],flankingvarnum)
	Vektoren=Data2Vectors(ExtSignatures)

		
# Postprocessing to weed out alternating variations
FinalSigs=[]

finalhisto=[0]*20
firsthisto=[0]*20

count=0

for z in range(len(ExtSignatures)):
	if CorrectedSigListen[z][-2]!=CorrectedSigListen[z][-1]:
		newsig=''
		for v in range(siglength):
			if CorrectedSigListen[z][-2][v]!=CorrectedSigListen[z][-1][v]:
				newsig+=Signatures[z][v]
			else:
				newsig+=CorrectedSigListen[z][-1][v]
		FinalSigs.append(newsig)
	else:
		FinalSigs.append(CorrectedSigListen[z][-1])			


# Writing out results:
string=''
for z in range(len(ExtSignatures)):
	string+=FinalSigs[z]+'\n'

f = open(SigPath+"_corrected",'w')
f.write(string)
f.close()

string=''
for z in range(len(ExtSignatures)):
	string+=CorrectedSigListen[z][0]+'\n'

f = open(SigPath+"_firstpass",'w')
f.write(string)
f.close()
				
# print "finalhisto:",					
# print finalhisto

# print "firsthisto:",
# print firsthisto

# print "unconverged:",
# print count

sys.exit()
