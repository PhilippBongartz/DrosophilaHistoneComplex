This repository contains the scripts, data and tools used to correct, assemble and test the drosophila histone complex. 

Some of these tools exist in several versions adapted to the different data sets. 

The file ‘Pipeline’ contains clean versions of all programs, implementing a pipeline capable of assembling a simulated complex, including a readme explaining how to run the pipeline and what each step is doing. 

The file ‘HistoneComplex’ contains those scripts specifically adapted to the histone data, as well as intermediate results such as the multiple sequence alignment, the trained neural networks, the clusters of signatures and the results of the graph touring. It also contains the final polished assembly sequence. 

The file ‘Transposons’ contains scripts adapted to the transposon data, as well as the transposon multiple sequence alignments, neural networks, signatures and ground truths of as many transposon data sets as fit into the repository. Additionally, it contains the data and scripts for the simulated data set used to benchmark the correction heuristic. 


README_Pipeline in /Pipeline describes all tools and scripts. In /Pipeline all data is created from scratch, starting with reads sampled from a simulated complex. In the following we describe the data and scripts that go beyond the pipeline and belong to specific data sets. 


Specifically adapted scripts:

TheanoNN.py in /HistoneComplex corrects the histone signatures.
Clustering.py in /HistoneComplex clusters the extended and corrected histone signatures.
GraphTouring.py in /HistoneComplex tours the graph constituted by the signature clusters. 

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
GraphTouring.py loads SavedClustersA2 and begins to construct a layered graph drawing, starting from the 5’ end of the complex.
It prints the ongoing distribution of clusters among layers and a comparison with the ground truth. It pickles the result under “Layer”.

InDels.py is an analysis step to add indels of several bases to the signatures as unified features. 

Transposons and simulated:
TransNNCorrector.py corrects the transposon and simulated signatures.
TransNNImprover2.py overfits the networks further to the transposon or simulated data if the Corrector didn’t learn enough. 
TransSimulated.py requires some additional tools and data. It is included to show how the simulated dataset was created.
Assessment2.py calculates the stats for the corrected transposon.
SimAssessment.py calculates the stats for the simulated signatures.


Data:

Multiple sequence alignments: *_MMA
Trained neural networks: Nets*
Signatures: uncorrected, first pass corrected, corrected.
TransposonCopies_* contain the ground truth obtained by clustering the flanking sequences. 
histoneseq.fasta is the template of the histone coding sequence
histoneconsensus.fa is the polished consensus of the final assembly of the whole complex
NewReadSeqinfointo contains information about how the signatures and unique sequences are arranged in the original reads
