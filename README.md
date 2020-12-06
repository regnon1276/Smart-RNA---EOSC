# Smart-RNA---EOSC
Abstract
Motivation: Infectious diseases from novel viruses are becoming a major public health concern. Increased interaction between human and nature, increased contact between human and wild animals increase the risk of apparition of new virus and the risk of interspecies transmission.
In this work, we have focused on two key aspects related to RNA and proteins evolution and interaction. 
Virus mutation is the first aspect by identifying  virus sequences mutation, tracking these sequence mutations then forecasting potential mutation area in the RNA sequence. 
We have combined the RNA alignment sequence method with a clustering algorithm and neural network in order to develop a novel method to insure the tracking and detect new location of the mutation.
Virus-Host interactions is the second key aspect, to forecast potential proteins and genes which can interact with a Virus and therefore being able to detect potential gene impacted and potential species impacted by the Virus.
Current computational prediction methods for novel viruses are based only on protein sequences. We have implemented a method allowing to combine protein sequence analysis together with Disease phenotypes (i.e., symptoms).

Results: We developed a tool and an application which has two major functionalities. The tool can  predict and recognize the origin of a Virus RNA sequence and the tool can allow it to detect new mutation locations in the sequence. Also the tool is able to perform protein alignment and homogeneity calculation, a method which is associated with a Phenotype analysis inspired by Natural language recognition techniques. 


The tool is composed of a python script proposed in this github repository and of a FPGA board ARTY which control the embedded AI.
The bistream and software to programm the board is also provided in this github.

The report done for EOSC project is showing also how to use the application.
