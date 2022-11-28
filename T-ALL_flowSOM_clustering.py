## I wrote this script for one of the undergraduates I work with to use. 
## The purpose of it is to cluster flow cytometry data using a self-organizing map called flowSOM. 
## He ran it on a couple example data sets and analyzed the clustering with my seurat wrappers. 



## FlowSOM code comes from https://github.com/Hatchin/FlowSOM
import sys
sys.path.append("/stor/home/jfm2773/opt/FlowSOM/script")
# For running the flowsom script locally
    ## Have to do this b/c PYPI package isn't updated
from cluster import *
from flowsom import *

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import matplotlib.pyplot as plt

import re
import pandas as pd
import numpy as np
import FlowCytometryTools
from FlowCytometryTools import test_data_dir, test_data_file
from FlowCytometryTools import FCMeasurement
from sklearn.cluster import AgglomerativeClustering
import os 
## This script needs to be run in the flowSOM virtual environment.



def flowSOM_wrapper(path, file_name, min_n, max_n, iter_n, prop_n, out_path, drop_cols= None):
    if drop_cols is None:
        drop= False
    tt = flowsom(os.path.join(path, file_name), if_fcs= False, if_drop= drop, drop_col= drop_cols)

    tt.som_mapping(50, # x_n: e.g. 100, the dimension of expected map
                   50, # y_n: e.g. 100, the dimension of expected map
                   tt.df.shape[1], # num of input columns
                   2.5, # sigma: e.g 2.5, the sigma of initialized weights
                   0.1, # lr: e.g 0.1, learning rate
                   500, # batch_size: 1000, iteration times
                   neighborhood= 'gaussian', # the initialized weights' distribution
                   tf_str=None, # Normalize beforehand
                   if_fcs=False, 
                   seed= 12345 ## set seed
                   ) 
    sample_name= re.sub(".csv", "", file_name)
    # som_output_weights = tt.map_som # SOM outputs
    # som_distance_map   = tt.map_som.distance_map() # distance map of SOM

    tt.meta_clustering(AgglomerativeClustering, #  cluster_class: e.g. KMeans, a cluster class, like "from sklearn.cluster import KMeans"
                       min_n, # the min proposed number of clusters
                       max_n, # the max proposed number of clusters
                       iter_n, # the iteration times for each number of clusters
                       resample_proportion=prop_n,
                       verbose=True 
                       )
    print("The best number of clusters is ", tt.bestk)
    tt.vis(t= 4, # the number of total nodes = t * bestk
           edge_color='b', 
           node_size=300, 
           with_labels=False)
    plt.savefig(os.path.join(out_path, "".join([sample_name, "_mst.png"])))

    tt.labeling()

    output_df = tt.df
    output_df.to_csv(os.path.join(out_path, "".join(['flowSOM_clustering_',file_name])), sep= ",")

    output_tf_df = tt.tf_df 
    output_tf_df.to_csv(os.path.join(out_path, "".join(['flowSOM_clustering_tf_', file_name])), sep= ",")

if __name__ == "__main__":
    input_path= "/stor/work/Ehrlich/T-ALL/myeloid/data/MyeloidComposition_Myeloid/csvs_from_fcs/normalized_data"
    out_path= "/stor/work/Ehrlich/T-ALL/myeloid/data/MyeloidComposition_Myeloid/clustering_output/"
    
    ##---------------------
    # Getting the csv files
    ##---------------------
    file_list= os.listdir(input_path)
    keyword= "Mouse1"
    file_list= [file for file in file_list if bool(re.search(keyword, file))]

    ##----------
    # Clustering
    ##----------
    ## the files that are passed to FlowSOM should be normalized and should only have the relevant markers. 
    ## If additional markers are present, specify them as a list to drop_cols argument.
    for file in file_list:
        flowSOM_wrapper(path      = input_path,
                        file_name = file,
                        out_path  = out_path,
                        drop_cols = None,
                        min_n     = 10, ## Min. num. of clusters to test
                        max_n     = 20, ## Max. num. of clusters to test
                        iter_n    = 1,  ## Num. of repeats for each cluster size
                        prop_n    = 0.6 ## Proportion of data that is resampled for each repeat
                        )
    ## Every time you change the clustering and want to save the output, make a new output directory.

    ## In the future, we might want to group the csv's together before clustering. 
    ## This will probably best be done as part of the pre-processing rather than in this script. 

   
