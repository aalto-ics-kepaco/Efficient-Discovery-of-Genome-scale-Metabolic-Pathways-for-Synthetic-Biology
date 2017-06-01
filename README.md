# Metabolic-Pathways-Identification


Implementation of the computational method kOpt, corresponding to the manuscript  "Efficient Discovery of Genome-scale Metabolic Pathways for Synthetic Biology", that has been submitted to Bionformatics journal.


**Abstract**

Motivation: Genome-scale metabolic networks and metabolic networks extended by new-to-nature reactions exceed the practical size of enumerating elementary flux modes, yet such enumeration is required to take advantage of the networks.
Results: A novel method, called kOpt, was developed and implemented to compute the first k flux solutions having increased flux pattern variability for reactions close to the product in the stoichiometric reaction network. The novel method proved feasible for genome-scale metabolic networks and compared favourably to a method enumerating the first k elementary flux modes. More specifically, the proposed method was approximately two orders of magnitude faster and yielded a more varied set of solutions.



**Contents**
-------------

* README.md 

  This file
 
* kOPt 

  This folder contains the MATLAB implementation of the kOPt method. Two versions of the code ara available, with the difference that one uses [cplex](http://www-03.ibm.com/software/products/en/ibmilogcpleoptistud) for the optimization process and the second one can be used with [GLPK](https://www.gnu.org/software/glpk/) optimization software. 
 
* kOptResults

  The results for each of the networks is stored in a separate mat file. The structure contains the optimization Objective values, running times, flux solutions and binary values of the flux solutions for optimization round. 

* EFMs

  The results for the enumeration of EFMs method. The running times and obtained EFMs for each network are saved in separate files.  

* EFMsResults

* data

 The data folder stores the stoichiometric (S) matrices representing the metabolic networks that we have used in the project. The S matrices are organized such that the rows represent metabolites, and the columns correspond to reactions in the network. 


* figures
  
  Matlab code for generation of the figures, and figures in different formats.

* calculateVariation 

  Scripts for calculation of the variability of the obtained solutions.

* supplementaryScripts
 
  Additional scripts required either for the methods or for the figure generation.

