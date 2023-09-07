# Weghted Ensemble method of Spatial Transcriptomics (WEST) 

## 1. Introduction 
WEST is a method that utilizes ensemble techniques to improve the performance and stability of spatial transcriptomics data analysis. 
It represents a significant advance in clustering spatial transcriptomics data, offering improved accuracy and flexibility compared to existing methods, 
making it a valuable tool for spatial transcriptomics data analysis.

## 2. Environment 
Python=3.8 
### Required packages: 

## 3.Tutorial 
The [tutorial](https://github.com/JiazhangCai/WEST/blob/main/tutorial.ipynb) provides a pipeline of how to implement
WEST on an actual sample, including the description for every parameter. 

## 4. Data used in the paper 
The DLPFC data used in the paper can be found [here](http://research.libd.org/spatialLIBD/). \n
The mouse embryo data used in the paper can be found [here](https://crukci.shinyapps.io/SpatialMouseAtlas/). \\ 
The simulated data used in the paper can be found in the folder [sim_data](https://github.com/JiazhangCai/WEST/tree/main/sim_data). \\ 
The reference-based synthetic data is already in the folder. \\  
The reference-free synthetic data can be generated from [shiny.rdata](https://github.com/JiazhangCai/WEST/blob/main/sim_data/ref_free/shiny.rdata) using 
the code [simulation_generate.R](https://github.com/JiazhangCai/WEST/blob/main/sim_data/ref_free/shiny.rdata).
