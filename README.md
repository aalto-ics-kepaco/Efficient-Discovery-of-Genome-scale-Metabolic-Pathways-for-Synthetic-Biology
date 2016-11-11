# Metabolic-Pathways-Identification


Master's thesis project on the topic "Identification of metabolic fluxes leading to the production of industrially relevant products". 

In this work we have developed and implemented a computational method, called kOpt, for enumeration of the k-first metabolic 
pathways with several biological constraints. Out methods is feasible for genome-scale metabolic networks as opposed to existing methods for enumeration of metabolic fluxes such as elementary flux modes and extreme pathways. We have compared the performance of kOpt against the method for enumeration of the first k elementary flux modes as developbed by Pey et al., in 2014, and the results demonstrate that our method is two orders of magnitude faster than the competing one. 

**Contents**
-------------

* README.md 

  This file
 
* kOPt 

  This folder contains the MATLAB implementation of the kOPt method. Two versions of the code ara available, with the difference that one uses [cplex](http://www-03.ibm.com/software/products/en/ibmilogcpleoptistud) for the optimization process and the second one can be used with [GLPK](https://www.gnu.org/software/glpk/) optimization software. 
 
* kOptResults

  The results

* EFMs

* EFMsResults

* data

 The data folder stores the stoichiometric (S) matrices representing the metabolic networks that we have used in the project. The S matrices are organized such that the rows represent metabolites, and the columns correspond to reactions in the network. 


* figures

* calculateVariation 

* supplementaryScripts

