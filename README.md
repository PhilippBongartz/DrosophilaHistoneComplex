# DrosophilaHistoneComplex
Data and code for the assembly of the drosophila histone complex

These are the scripts, data and tools used to correct, assemble and test the histone complex. 


Tools:

PW_ReAligner.c realigns sequences to a multiple sequence alignment, progressively refining the alignment.
Correlation.c calculates the statistical significance of the correlation between columns in a multiple sequence alignment. 

Both are run by piping the multiple sequence alignment MSA into the tool: cat MSA | ./tool

PW_ReAligner has no dependencies. Correlation requires the gnu scientific library:

gcc -Wall -I/usr/local/include -c Correlation.c
gcc -L/usr/local/lib Correlation.o -lgsl -lgslcblas -lm -o Correlation


Scripts:

TheanoNN.py corrects the histone signatures.
Clustering.py clusters the extended signatures.
GraphTouring.py tours the graph constituted by the signature clusters. 

These three are the essential steps for the Histone complex assembly, so I’ll describe them in some more detail:

python TheanoNN.py UltimateVariations2 
TheanoNN.py takes "UltimateVariations2" as an argument and creates the "FirstUltimateVariations" and "FinalUltimateVariations"
The neural networks are supposed to be overfitted to the data. For some nets this requires several rounds of training. 
If the initial learning rate is high, the nets learn fast, but some might become stuck in local minima. 
If a net seems to be stuck, it is best to just delete it and let the next round of learning initialize a new net.

python Clustering.py 
Clustering.py loads the corrected Signatures and additional information about the reads and clusters signature triplets.
These Clusters are pickled under SavedClustersA2.

python GraphTouring.py 
GraphTouring.py loads SavedClustersA2 and begins to construct a layered graph drawing, starting from the 5’ end of the complex. It prints the ongoing distribution of clusters among layers and a comparison with the ground truth. It pickles the result under “Layer”.

InDels.py is an analysis step to add indels of several bases to the signatures as unified features. 

Transposons and simulated:
TransNNCorrector.py corrects the transposon and simulated signatures.
TransNNImprover2.py overfits the networks further to the transposon or simulated data if the Corrector didn’t learn enough. 
TransSimulated.py requires some additional tools and data. It is included to show how the simulated dataset was created.
Assessment2.py calculates the stats for the corrected transposon.
SimAssessment.py calculates the stats for the simulated signatures.


Data:

Multiple sequence alignments: *MMA
Trained neural networks: Nets*
Signatures: uncorrected, first pass corrected, corrected.
TransposonCopies* contain the ground truth obtained by clustering the flanking sequences. 

