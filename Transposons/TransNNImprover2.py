#!/usr/bin/env python
# coding: latin1

import random
import numpy as np
import time
import sys
import os
import operator
import copy
import pickle
import cPickle
import theano
import theano.tensor as T
from theano.tensor.nnet import conv
from theano.tensor.nnet import softmax
from theano.tensor import shared_randomstreams
from theano.tensor.signal import downsample

### Der NN-Corrector für die Transposon-Signaturen
# 
# Diese Version läd bereits trainierte Netze und lässt die schlechten nochmal mehr Lernen. 
# 
### 

# Code für neural network stuff
def linear(z): return z
def ReLU(z): return T.maximum(0.0, z)
from theano.tensor.nnet import sigmoid
from theano.tensor import tanh
from theano.tensor.nnet import softmax as rowsoftmax

def softmax(x):
	return rowsoftmax(x.T).T

rng = np.random.RandomState(1234)
srng = T.shared_randomstreams.RandomStreams(rng.randint(999999))

def drop(inputv, p=0.5, rng=rng): 
	"""
	:type input: numpy.array
	:param input: layer or weight matrix on which dropout resp. dropconnect is applied
	
	:type p: float or double between 0. and 1. 
	:param p: p probability of NOT dropping out a unit or connection, therefore (1.-p) is the drop rate.
	
	"""            
	mask = srng.binomial(n=1, p=p, size=inputv.shape, dtype=theano.config.floatX)
	return inputv * mask

def norm(inputv, p=0.5):
	return inputv * p


class Network(object):

	def __init__(self,layersizes, activation_fn1=sigmoid, activation_fn2=softmax):
		if len(layersizes)!=3:
			print "Dies ist ein DreierNN."
			sys.exit()

		self.sizes=layersizes
		self.activation_fn1=activation_fn1
		self.activation_fn2=activation_fn2
		self.WeightInitializer()

		# Input, Target und Output sind jeweils Vektoren.
		self.x = T.dmatrix("x")  
		self.y = T.dmatrix("y")
		
		self.a1 = T.dmatrix("a1") 
		self.a2 = T.dmatrix("a2") 

		self.z0 = T.dmatrix("z0") 
		self.z1 = T.dmatrix("z1")   # Das ist output


	def WeightInitializer(self):	
		self.W0 = theano.shared(
			np.asarray(
				np.random.normal(loc=0.0, scale=np.sqrt(5.0/self.sizes[0]), size=(self.sizes[1],self.sizes[0])),
					dtype=theano.config.floatX),
			name='W0', borrow=True)

		self.W1 = theano.shared(
			np.asarray(
				np.random.normal(loc=0.0, scale=np.sqrt(1.0/self.sizes[1]), size=(self.sizes[2],self.sizes[1])),
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


	def SGD(self, minibatchsize, trainingsdata, epochs, eta, lambd, dropout=0):
		# Hier wird die Costfunktion in Abhängigkeit von self.x, self.y und den Weights und Biases berechnet.
		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		# Mit Dropout:
		if dropout:
			self.z11 = drop(self.z1)*2.0  # Das ist z1 mit gedroppten Einträgen.
			self.a2 = T.dot(self.W1, self.z11) + self.B1

		# Ohne Dropout:
		if not dropout:
			self.a2 = T.dot(self.W1, self.z1) + self.B1

		self.z2 = self.activation_fn2(self.a2)

		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd
		qcost = (((self.z2 - self.y)**2).sum()) * 0.5 + regularization 

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

		for i in range(epochs):
			z = random.randint(0,len(trainingsdata)-minibatchsize)

			mini_batch = trainingsdata[z:z+minibatchsize]

			changes = [0.0 for param in self.params] # Upgedated wird am Ende des Minibatch
			for (x,y) in mini_batch:
				upgrads = gradients(x,y)
				changes = [ change - upgrad*eta/float(minibatchsize) for upgrad,change in zip(upgrads,changes)]   # /minibatchsize damit die Lernrate davon unabhängig bleibt.

			#self.params = [ param-change for param, change in zip(self.params, changes) ]
			updaterW0(changes[0])
			updaterW1(changes[1])
			updaterB0(changes[2])
			updaterB1(changes[3])


	def Accuracy(self,trainingsdata, lambd):		
		# Hier wird die Costfunktion in Abhängigkeit von self.x, self.y und den Weights und Biases berechnet.
		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)
		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)
		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd
		qcost = (((self.z2 - self.y)**2).sum()) * 0.5 + regularization 
		#qcost = -T.mean(self.y*T.log(self.z2)+(1-self.y)*T.log(1-self.z2)) + regularization

		starget = T.argmax(self.y,axis=0)
		y_out = T.argmax(self.z2, axis=0)
		grads = T.grad(qcost, self.params)

		cost = theano.function([self.x,self.y],qcost)
		feedforward = theano.function([self.x],y_out)
		gradients = theano.function([self.x,self.y],grads)
		target = theano.function([self.y],starget)

		output = theano.function([self.x],self.z2)

		# correctpositives, falsepositives, correctnegatives, falsenegatives und cost
		avcost=0.0
		correct=0.0
		count=0.0

		Targets=[0,0,0,0,0]
		Predicts=[0,0,0,0,0]

		for (x,y) in trainingsdata:
			avcost+=cost(x,y)

			ff=feedforward(x)
			tt=target(y)

			Targets[tt]+=1
			Predicts[ff]+=1

			if ff==tt:
				correct+=1.0
	
		print "Average Cost = {}.".format(avcost/float(len(trainingsdata)))
		print "Correct: {} prozent.".format(float(int(1000.0*correct/float(len(trainingsdata))))/10.0)	
		print "{} -> {}".format(Targets,Predicts)
		return (avcost/float(len(trainingsdata)),float(int(1000.0*correct/float(len(trainingsdata))))/10.0,Targets,Predicts)


	def PredictedBase(self, x):		

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		y_out = T.argmax(self.z2, axis=0)
		
		feedforward = theano.function([self.x],y_out)

		y=feedforward(x)
		Index2Basen=['a','c','g','t','-']

		return Index2Basen[y]

	def PredictedBases(self, InputData):		
		# Hier wird die Costfunktion in Abhängigkeit von self.x, self.y und den Weights und Biases berechnet.

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		y_out = T.argmax(self.z2, axis=0)
		
		feedforward = theano.function([self.x],y_out)
		Index2Basen=['a','c','g','t','-']
		NewBases=[]

		for x in InputData:
			y=feedforward(x)
			NewBases.append(Index2Basen[y])

		return NewBases
		
# Code für Data stuff

# Hier mal Vectoren ohne Auslassung: Die Targets werden beim Training dann ignoriert.
Basen2Index={'a':0,'c':1,'g':2,'t':3,'-':4}
def Data2Vectors(ISignatures):	
	Sigvectoren=[]
	for i in range(len(ISignatures)):
		OneInput=np.zeros((len(ISignatures[0])*5, 1))   # Zweiervector pro Variation
		
		for v in range(len(ISignatures[0])):
			if ISignatures[i][v]!=' ':
				OneInput[Basen2Index[ISignatures[i][v]]+v*5]=1.0
		Sigvectoren.append(OneInput)	
	return Sigvectoren	




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

	print VarPATH

	Signatures=[]
	f=open(VarPATH,'r')

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
	print "siglength: {}".format(siglength)

### MAIN ###

eta=1.01
lambd=1.0 / float((siglength-1)*2)**2

Vektoren=Data2Vectors(Signatures)
batchsize=500

nichtsgelernt=0

Nets=[]
for Var in range(0,siglength):
	namestring="Var_"+str(Var)+"_net"
		
	#if 1:
	if numero>-1 and os.path.exists("Nets_"+str(numero)+"/"+namestring):
		trainingsdata=[(np.concatenate((Vektor[:Var*5],Vektor[(Var+1)*5:]),axis=0),Vektor[Var*5:(Var+1)*5]) for Vektor in Vektoren if np.linalg.norm(Vektor[Var*5:(Var+1)*5])>0.5]   # Var
		batchsize=min(batchsize,len(trainingsdata))
		net = cPickle.load(open("Nets_"+str(numero)+"/"+namestring,"rb"))		
		# net = Network([(len(Signatures[0])-1)*5,10,5])
		print "Var: {}, {}, bisher nichts gelernt {}".format(Var,VarPATH, nichtsgelernt)
		z=0
		Predicts=[0,0,0,0,0]
		weiter=True
		lastcorrect=0
		lastcost=1.0
		start=time.time()
		(avcost,correct,Targets,Predicts)=net.Accuracy(trainingsdata, lambd)

		# Das hier füge ich für diese Version ein: Jedes Netz wird ein bisschen overfitted:
		net.SGD(50,trainingsdata,5000,eta,lambd,0)
		(avcost,correct,Targets,Predicts)=net.Accuracy(trainingsdata, lambd)

		maxavcost=avcost
		eta=0.3
		if sum([1 for x in range(5) if Predicts[x]>20])<2:  # sum(Predicts)==max(Predicts):   
			while (sum([1 for x in range(5) if Predicts[x]>20])<2 or avcost>0.05) and z<10:
				net.SGD(50,trainingsdata,5000,eta,lambd,0)
				(avcost,correct,Targets,Predicts)=net.Accuracy(trainingsdata, lambd)
				if avcost<maxavcost:
					maxavcost=avcost
					maxNet=copy.deepcopy(net)			
				if lastcost<avcost+0.001:
					eta/=2.0
				z+=1
				lastcorrect=correct	
				lastcost=avcost

		if sum([1 for x in range(5) if Predicts[x]>20])<2:
			nichtsgelernt+=1		

		maxNet=copy.deepcopy(net)
		maxavcost=avcost		

		print "{} sec trainiert.".format(int(time.time()-start))
		print
	
		cPickle.dump(maxNet,open("Nets_"+str(numero)+"/"+namestring,"w"))

	if numero==-1 and os.path.exists("Nets_Sim/"+namestring):
		trainingsdata=[(np.concatenate((Vektor[:Var*5],Vektor[(Var+1)*5:]),axis=0),Vektor[Var*5:(Var+1)*5]) for Vektor in Vektoren if np.linalg.norm(Vektor[Var*5:(Var+1)*5])>0.5]   # Var
		# net = Network([(len(Signatures[0])-1)*5,10,5])
		net = cPickle.load(open("Nets_Sim/"+namestring,"rb"))
		print "Var: {}".format(Var)
		z=0
		Predicts=[0,0,0,0,0]
		weiter=True
		lastcorrect=0
		lastcost=1.0
		start=time.time()
		net.SGD(batchsize,trainingsdata,1000,10.0,lambd,1)
		(avcost,correct,Targets,Predicts)=net.Accuracy(trainingsdata, lambd)
		eta=1.0
		while sum(Predicts)!=max(Predicts) and z<5:
			net.SGD(50,trainingsdata,5000,eta,lambd,0)
			(avcost,correct,Targets,Predicts)=net.Accuracy(trainingsdata, lambd)
			if avcost<maxavcost:
				maxavcost=avcost
				maxNet=copy.deepcopy(net)			
			if lastcost<avcost+0.001:
				eta/=2.0
			z+=1
			lastcorrect=correct	
			lastcost=avcost

		maxNet=copy.deepcopy(net)
		maxavcost=avcost		

		print "{} sec trainiert.".format(int(time.time()-start))
		print
	
		cPickle.dump(maxNet,open("Nets_Sim/"+namestring,"w"))			


sys.exit()




