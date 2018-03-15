# Simulating-and-visualizing-Boolean-models

## Introduction

These are scripts to create Boolean models from prior knowledge interaction networks and network optimization based on experimental data stored in the standardized input format of [MIDAS tables](https://academic.oup.com/bioinformatics/article/24/6/840/193440).
The network optimization tool that was used here is the R package [CellNOptR](http://www.cellnopt.org). The logical models obtained from network optimization are further tested by the BoolNet R package.

...material in process of uploading

#########

## File description

The files  
*kras_model_conceptual.cfg*  
*kras_model_conceptual.bnd*  
are model files of a protein interaction network for the stochastic Boolean simulator platform [MaBoSS](http://maboss.curie.fr).

The R script  
*maboss_wrapper_plotting.R*  
is a pipeline of functions to introduce mutations and inhibitions in this model (or any other in the MaBoSS format), perform simulations and subsequently to load and visualize the simulation results.
The functions themselves are contained in *maboss_wrapper_plotting_functions.R*.
