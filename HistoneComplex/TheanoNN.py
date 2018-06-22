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
# It takes "UltimateVariations2" as an argument and creates the "FirstUltimateVariations" and "FinalUltimateVariations"
# The neural networks are supposed to be overfitted to the data. For some nets this requires several rounds of training. 
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

# Reading in information about the placement of signatures and large-scale vars in reads

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

# This function extends signatures by adding bases of the neighbouring signatures:
def RealDataSignatureExtension(Signatures, flankingvarnum):
	ExtSignatures=[]
	KlassPattern={'Li':"acgt",'Re':"tgca",'Ze':"aaaa",'DA':"cccc",'AN':"gggg",'NN':"tttt"}
	for x in range(len(IndexReads)):
		for y in range(len(IndexReads[x])):
			if KlassReads[x][y] in ['No','Ku']:	

				sig=Signatures[IndexReads[x][y]]
				linkssig=''
				rechtssig=''

				for z in range(y+1,len(IndexReads[x])):
					if KlassReads[x][z] in ['No','Ku']:
						rechtssig+=Signatures[IndexReads[x][z]]
					elif KlassReads[x][z] in KlassPattern:
						klasspatt=KlassPattern[KlassReads[x][z]]*siglength
						rechtssig+=klasspatt[:siglength]					
					else:
						rechtssig+=' '*siglength
				rechtssig+=' '*flankingvarnum		

				for z in range(y-1,-1,-1):
					if KlassReads[x][z] in ['No','Ku']:
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

		
# Reading in the signatures 
if len(sys.argv)>1 and os.path.exists(sys.argv[1]):
	VarPATH=sys.argv[1]
	print VarPATH

	Signatures=[]
	f=open(VarPATH,'r')

	for line in f.readlines():
		if line[-1]=='\n':
			line=line[:len(line)-1]

		if line == ' '*len(line):	
			Signatures.append(' ')
		else:	
			Signatures.append(line)

	for sig in Signatures:
		if len(sig)>1:
			siglength=len(sig)	

	for x in range(len(Signatures)):
		if Signatures[x]==' ':
			Signatures[x]=siglength*' '
	count=0
	for x in range(len(Signatures)-1):
		if len(Signatures[x])!=siglength:
			count+=1
		if len(Signatures[x])!=len(Signatures[x+1]):
			print len(Signatures[x+1])
			print Signatures[x+1]

	print "count",
	print count		

	for sig in Signatures:
		if len(sig)!=siglength:
			print "Siglengtherror"
			sys.exit()		

	Variations=range(siglength)
	n=siglength-1

	print "siglength: {}".format(siglength)

	#This parameter determines how many bases of the neighbouring signatures are used for the correction.
	flankingvarnum=siglength

	length=siglength+flankingvarnum+flankingvarnum

	ExtSignatures=RealDataSignatureExtension(Signatures,flankingvarnum)
	
	#Truedata=TrueData(Signatures)

	for x in range(len(ExtSignatures)-1):
		if len(ExtSignatures[x])!=len(ExtSignatures[x+1]):
			print len(ExtSignatures[x+1])


	print "length:",
	print length
	print len(ExtSignatures[0])

else: 
	print "Please provide a signature file."


#########################################################
			
# Signatures to vectors
def Data2Vectors(ISignatures):	
	Sigvectoren=[]
	for i in range(len(ISignatures)):
		OneInput=np.zeros((len(ISignatures[0])*5, 1))   # five entries per variation
		
		for v in range(len(ISignatures[0])):
			if ISignatures[i][v]!=' ':
				OneInput[Basen2Index[ISignatures[i][v]]+v*5]=1.0
		Sigvectoren.append(OneInput)	
	return Sigvectoren		


# This is the neural network implementation: 
class Network(object):

	def __init__(self,layersizes, activation_fn1=sigmoid, activation_fn2=sigmoid):
		if len(layersizes)!=3:
			print "Dies ist ein DreierNN."
			sys.exit()

		self.sizes=layersizes
		self.activation_fn1=activation_fn1
		self.activation_fn2=activation_fn2
		self.WeightInitializer()

		# Input and target
		self.x = T.dmatrix("x")  
		self.y = T.dmatrix("y")
		
		self.a1 = T.dmatrix("a1") 
		self.a2 = T.dmatrix("a2") 

		self.z0 = T.dmatrix("z0") 
		self.z1 = T.dmatrix("z1")   # the output


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

		# The neural network computation graph
		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		# The regularization
		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd

		act_reg = ((self.z1-0.5)**2).sum()

		# quadratic cost + regularisation
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

		# The learning: 
		for i in range(epochs):
			z = random.randint(0,len(trainingsdata)-minibatchsize)

			mini_batch = trainingsdata[z:z+minibatchsize]

			changes = [0.0 for param in self.params] # Updating at the end of a minibatch
			for (x,y) in mini_batch:
				upgrads = gradients(x,y)
				changes = [ change - upgrad*eta/float(minibatchsize) for upgrad,change in zip(upgrads,changes)]   # /minibatchsize damit die Lernrate davon unabhängig bleibt.

			#self.params = [ param-change for param, change in zip(self.params, changes) ]
			updaterW0(changes[0])
			updaterW1(changes[1])
			updaterB0(changes[2])
			updaterB1(changes[3])

	# This is experimental code which adapts the learning rate for faster conversion:
	def SelftuningSGD(self, minibatchsize, trainingsdata, epochs, eta, lambd, selftuningstep):		

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

			changes = [0.0 for param in self.params] # Upgedating at the end of a minibatch
			for (x,y) in mini_batch:
				upgrads = gradients(x,y)
				changes = [ change - upgrad*eta/float(minibatchsize) for upgrad,change in zip(upgrads,changes)]  

			#self.params = [ param-change for param, change in zip(self.params, changes) ]
			updaterW0(changes[0])
			updaterW1(changes[1])
			updaterB0(changes[2])
			updaterB1(changes[3])

			# tuning eta via false positive and true positives
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
						if not trainingsdatacorrect[yyy]:   # not yet correctly corrected
							neuverbessert+=1
						trainingsdatacorrect[yyy]=1
					else:		# not correctly corrected
						if trainingsdatacorrect[yyy]:  # had been correctly corrected
							neuverschlechtert+=1
						trainingsdatacorrect[yyy]=0	

				if neuverschlechtert==0:
					eta*=1.1
				elif 2*neuverschlechtert>neuverbessert:
					eta/=1.1	

				#print "Average Cost = {}.".format(avcost/float(len(trainingsdata)))
				print "Correct: {} prozent. {} neuverbessert, {} neuverschlechtert. Neu-Eta: {}".format(float(int(1000.0*correct/float(len(trainingsdata))))/10.0, neuverbessert, neuverschlechtert,eta)					
				self.Accuracy(trainingsdata,lambd)



	# Output of the current accuracy and error
	def Accuracy(self,trainingsdata, lambd):		

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd

		#qcost = (((self.z2 - self.y)**2).sum())**0.5 + regularization #+ act_reg

		# quadratic cost
		qcost = (((self.z2 - self.y)**2).sum()) * 0.5 + regularization #+ act_reg
		
		starget = T.argmax(self.y,axis=0)
		y_out = T.argmax(self.z2, axis=0)
		grads = T.grad(qcost, self.params)

		cost = theano.function([self.x,self.y],qcost)
		feedforward = theano.function([self.x],y_out)
		gradients = theano.function([self.x,self.y],grads)
		target = theano.function([self.y],starget)

		output = theano.function([self.x],self.z2)

		# hiddenstate = theano.function([self.x],self.z1)

		avcost=0.0
		correct=0.0

		basehisto=[0]*self.sizes[2]
		truehisto=[0]*self.sizes[2]
		truehits=[0]*self.sizes[2]

		for (x,y) in trainingsdata:
			avcost+=cost(x,y)

			ff=feedforward(x)
			tt=target(y)

			basehisto[ff]+=1
			truehisto[tt]+=1

			if ff==tt:
				truehits[tt]+=1
				correct+=1.0

		#print "pred. {} vs truth {}: truehits {}".format(basehisto,truehisto,truehits)
		print "Average Cost = {}.".format(avcost/float(len(trainingsdata)))
		print "Correct: {} prozent.".format(float(int(1000.0*correct/float(len(trainingsdata))))/10.0)	
		return float(int(1000.0*correct/float(len(trainingsdata))))/10.0
		#return avcost/float(len(trainingsdata))

	# Predicting one base
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

	# Predicting many bases
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
			NewBases.append(Index2Basen[y])
			

		return NewBases
		
	# A diagnostic function to detect saturated neurons.
	def AvAct(self,trainingsdata, lambd):		

		self.a1 = T.dot(self.W0, self.x) + self.B0
		self.z1 = self.activation_fn1(self.a1)

		self.a2 = T.dot(self.W1, self.z1) + self.B1
		self.z2 = self.activation_fn2(self.a2)

		regularization = ((self.W0**2).sum() + (self.W1**2).sum() + (self.B0**2).sum() + (self.B1**2).sum())*lambd

		#qcost = (((self.z2 - self.y)**2).sum())**0.5 + regularization #+ act_reg

		# quadratic cost
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

		# hiddenstate = theano.function([self.x],self.z1)

		avcost=0.0
		correct=0.0

		basehisto=[0,0,0,0,0]

		# av. activation to detect neuron saturation
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



################ Calculating and pickling all nets

for v in range(siglength):  #(125,siglength): #,siglength):
	print
	print
	print "Variation {}:".format(v)
	Var=v+flankingvarnum
	trainingsdata=[(np.concatenate((Vektor[:Var*5],Vektor[(Var+1)*5:]),axis=0),Vektor[Var*5:(Var+1)*5]) for Vektor in Vektoren if np.linalg.norm(Vektor[Var*5:(Var+1)*5])>0.5]   # Var

	namestring="UltVar"+str(v)+"fvn"+str(flankingvarnum)+"_net"

	if os.path.exists("Nets2/"+namestring):
		net=cPickle.load(open("Nets2/"+namestring,"rb"))
	else:
		net = Network([(length-1)*5,5,5], tanh, sigmoid)

	lambd=1.0 / float((length-1)*5)

	acc=net.Accuracy(trainingsdata, lambd)

	if acc < 95:
		print "Further Training:"
		eta=0.0030000
		net.SGD(10,trainingsdata,500000,eta, lambd)
		net.Accuracy(trainingsdata, lambd)		


	namestring="UltVar"+str(v)+"fvn"+str(flankingvarnum)+"_net"
	cPickle.dump(net,open("Nets2/"+namestring,"w"))

#sys.exit()

################# Correction:
# Loading the nets and creating new signatures
CorrectedSigListen=[[] for z in range(len(Signatures))] # several runs of correction

print "Loading nets."
Nets=[]
for v in range(siglength):  #(125,siglength): #,siglength):
	namestring="UltVar"+str(v)+"fvn"+str(flankingvarnum)+"_net"
	net=cPickle.load(open("Nets2/"+namestring,"rb"))
	Nets.append(net)

for d in range(5):  # Correction
	print "Round {}.".format(d)
	for z in range(len(Signatures)):
		CorrectedSigListen[z].append('')

	for v in range(siglength):
		print v,
		sys.stdout.flush()

		Var=v+flankingvarnum
		InputData=[ np.concatenate((Vektor[:Var*5],Vektor[(Var+1)*5:]),axis=0) for Vektor in Vektoren]
		
		NewBases=Nets[v].PredictedBases(InputData)

		for z in range(len(Signatures)):
		
			if Signatures[z][v]==' ':
				CorrectedSigListen[z][-1]+=' '
			else:
				CorrectedSigListen[z][-1]+=NewBases[z]
			
	print		
	# At the end of each run new vectors are calculated for newsigs:
	ExtSignatures=RealDataSignatureExtension([SigListe[d] for SigListe in CorrectedSigListen],flankingvarnum)
	Vektoren=Data2Vectors(ExtSignatures)

		
# Postprocessing to weed out alternating variations
FinalSigs=[]

finalhisto=[0]*20
firsthisto=[0]*20

count=0

for z in range(len(Signatures)):
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
for z in range(len(Signatures)):
	string+=FinalSigs[z]+'\n'

f = open('FinalUltimateSigs','w')
f.write(string)
f.close()

string=''
for z in range(len(Signatures)):
	string+=CorrectedSigListen[z][0]+'\n'

f = open('FirstUltimateSigs','w')
f.write(string)
f.close()
				
print "finalhisto:",					
print finalhisto

print "firsthisto:",
print firsthisto

print "unconverged:",
print count

sys.exit()
