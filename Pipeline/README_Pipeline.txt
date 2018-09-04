README_Pipeline

The pipeline consists of ten scripts and tools that implement the full assembly of a simulated repeat complex. For the assembly of a real complex, additional manual analysis steps are necessary to identify unique large scale variations.

All C programs have a -h flag, that explains usage. 

All python script use python 2.7.

The purpose of this pipeline is to illustrate the workflow of a repeat complex assembly as described in “Deep Repeat Resolution - The Assembly of the Drosophila Histone Complex”. But of course, the tools can also be used in isolation. Each tool has a default output path. It is possible to override this output path, but the last parts of the pipeline, that are quite specific to the simulated complex, have the default paths of input data hardcoded. So if the entire pipeline is run, it is better to just use the default paths. 

The whole pipeline has been tested on mac and linux. As development is ongoing and each pipeline step depends on its predecessors, there is nonetheless a lot of room for things to go wrong. If you experience problems running the pipeline, please do not hesitate to mail philipp.bongartz@h-its.org with a description of the problem. 

The pipeline can also be applied to real data, starting with a fasta file containing PacBio reads sampled from a repeat complex. However, ideally in a preprocessing step all unique sequences in the data would be detected, classified and excised, so that only the repeat sequences are arranged into the multiple sequence alignment. As mentioned in the paper, this step is currently being done manually and its automatisation is a possible direction for future work.

In the following I provide a description of every program in the pipeline, each closing with the compile and run commands. 



SimulatedComplex.py

Simulates a repeat complex with parameters close to the histone complex: repeat length 5kbp, copy number 100, SNP number 200 and coverage 90. It creates ’SimulatedTemplate.fasta’ the repeat template, ‘SimulatedReads.fasta’ reads sampled from the simulated complex and ‘ReadPlacements’ which contain the ground truth read placements within the complex. 

python SimulatedComplex.py



LargeScaleVars.c

Cuts reads into instances of the repeat sequence. Aditionally prints stats about the order and distance of sections of the repeat sequence, which is a first step to identify many large scale insertions, deletions and duplications. It creates ‘SimulatedSeq.fasta’ which contains the sequences into which the reads have been cut and ‘SimulatedReadSeqInfo’ which contains the information about which read has been cut into which sequences. 

gcc LargeScaleVars.c -o LargeScaleVars -lm
./LargeScaleVars SimulatedTemplate.fasta SimulatedReads.fasta



InitialAligner.c

Creates a multiple sequence alignment ‘SimulatedMSA’ from the sequences by aligning them to the repeat template. This is a very rough first muli-alignment, which still has to be refined. It also outputs ‘SimulatedSeqClass’. This file divides sequences into repeat sequences ‘r’, which are arranged into the initial MSA and deviating sequences ‘l’, which will be separately processed. 

gcc InitialAligner.c -o InitialAligner
./InitialAligner SimulatedTemplate.fasta SimulatedSeq.fasta



PW_ReAligner.c

Refines the initial multiple sequence alignment ‘SimulatedMSA’ by aligning rows iteratively to the multi-alignment until convergence of the sum of pairwise scores. Outputs ‘SimulatedMSAreal’.

gcc -mcmodel=medium PW_ReAligner.c -o PW_ReAligner
./PW_ReAligner SimulatedMSA



Correlation.c

Calculates the statistical significance of intersections between base groups in the columns of ‘SimulatedMSAreal’. Outputs ‘SimulatedSignatures’. Requires gnu scientific library. 

gcc -Wall -std=c99 -I/usr/local/include -c Correlation.c
gcc -L/usr/local/lib Correlation.o -lgsl -lgslcblas -lm -o Correlation
./Correlation SimulatedMSAreal



SeqClustering.c

Just to make the Pipeline complete, we implemented this simple version of a sequence clustering algorithm. It is used to cluster the sequences that did not fit the template in InitialAligner according to alignment score. For the simulated data set this means disambiguating the unique flanking sequences. The output ‘SeqClusters’ is then used in the graph touring to add those clusters to the clusters of repeat sequences calculated by Clustering.py

gcc -Wall SeqClustering.c -o SeqClustering
./SeqClustering SimulatedSeq.fasta SimulatedSeqClass



TheanoNN.py

Trains a neural network for every entry in the signatures ‘SimulatedSignatures’. Uses these networks to correct the signatures and outputs ‘SimulatedSignatures_corrected’ and ‘SimulatedSignatures_firstpass’. Can (and should) be run several times, nets are reloaded and trained further if the accuracy is still low. This can be done via the parameter <start_variation> that can be used to run TheanoNN.py on several cores going through all variations on all cores but starting with different variations to avoid collision. If the spacing x is chosen too small, collisions might crash the script. 

python TheanoNN.py SimulatedSignatures SimulatedSeqClass SimulatedReadSeqInfo 0
python TheanoNN.py SimulatedSignatures SimulatedSeqClass SimulatedReadSeqInfo 0+x
python TheanoNN.py SimulatedSignatures SimulatedSeqClass SimulatedReadSeqInfo 0+2x
…



Correction.py

This the part of TheanoNN.py that loads all nets and does the correction of Signatures. TheanoNN.py only corrects the signatures after it has gone through training all nets. If it is run on many cores and the start_variation is spaced appropriately all nets will have received some amount of training before each TheanoNN.py has gone through all of them. Correction.py can be used once all nets are created.

python Correction.py SimulatedSignatures SimulatedSeqClass SimulatedReadSeqInfo



Clustering.py

Clusters extended ‘SimulatedSignatures_corrected’ and outputs ‘Clusters’ and ‘Centroids’. The clustering is relatively simple and hopefully robust. The only internal ‘magic’ number is that we choose the centroid coverage such that more then 10% of the extended signatures are initial centroids. 

python Clustering.py SimulatedSignatures_corrected SimulatedSeqClass SimulatedReadSeqInfo



GraphTouring.py

Uses ‘Clusters’ and ‘Centroids’ to sort the signature clusters into layers that constitute copy groups of the assembly. Adds ‘SeqCluster’ to ideally recreate the simulated complex from flanking sequence to flanking sequence. Furthermore it loads 'SimulatedSignatures_corrected’, 'SimulatedSeqClass' and 'SimulatedReadSeqInfo' to create the extended signatures that are the basis of the ‘Clusters’ and ‘Centroids’. It also calculates the ground truth based on ‘ReadPlacements’, 'SimulatedReads.fasta’ and 'SimulatedSeq.fasta' to give a first assessment of the graph touring results. 

python GraphTouring.py



Results.py

Creates some stats about the correction and the assembly. It loads most of the data also used in GraphTouring.py and calculates error rate reductions. It also prints out how accurately the ground truth complex has been reconstructed by the clustering + graph touring. 

python Results.py



