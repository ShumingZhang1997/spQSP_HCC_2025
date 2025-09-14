#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 3 02:07:32 2024

@author: szhang
"""

from dask.distributed import Client, get_worker
from dask_jobqueue import SLURMCluster
import os
import numpy as np 

import matplotlib.pyplot as plt
from pyabc import ABCSMC, Distribution, RV, sampler
from pyabc import LocalTransition, MedianEpsilon
from pyabc.visualization import plot_data_callback, plot_kde_2d, plot_kde_1d

import pcdl
import shutil
from SPQSP_Model import SPQSP_Model

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import pandas as pd

############ PhysiCell Model ############
execFileName = "/home/szhan121/scratch4-apopel1/szhan121/MARCC/HCC/linux/HCC_s_sim"
model = SPQSP_Model(execFile = execFileName, time = 1, grid = 3, grid_interval= 2)

def model_summary(output_folders, metrics):
    #os.environ['R_HOME'] = "/data/apps/extern/spack_on/gcc/9.3.0/r/4.3.0-5w26bclf7ogtbbwwfx337icyt6p53t3m/rlib/R"
    # Initialize R to Python conversion
    pandas2ri.activate()
    abs_path = os.path.abspath(output_folders)
    os.chdir(abs_path)
    # Run the R script file
    R_script = '/home/szhan121/scratch4-apopel1/szhan121/MARCC/HCC/postprocessing/cell_kernel_marcc.R'
    robjects.r('.libPaths(c("/home/szhan121/R/x86_64-pc-linux-gnu-library/4.3", .libPaths()))')

    # Source the R script properly
    robjects.r(f'source("{R_script}")')

    # Now access R variables in Python
    try:
        HCC_response_data = pd.DataFrame(robjects.r['spatial_merged_df']).T  # Fetch the variable
        HCC_response_data_colname = list(robjects.r('colnames(spatial_merged_df)'))
        HCC_response_data.columns = HCC_response_data_colname
        # Debugging: Print the DataFrame and its columns
        print("DataFrame columns:", HCC_response_data.columns)
        print("DataFrame head:\n", HCC_response_data.head())
        spatial_metric = HCC_response_data[metrics].mean()
        return spatial_metric
    except:
        print(f"Variable '{metrics}' not found in R environment")
        # List available R objects
        print("Available R objects:")
        print(robjects.r('ls()'))
        return 0

def RunSPQSPmodel(pars):
    # Set the parameters in PhysiCell class
    SampleID = os.getpid()  # worker number
    Parameters = {"TCell_moveSteps": np.array([pars["TCell_moveSteps"]]),
                  "Fib_moveSteps": np.array([pars["Fib_moveSteps"]]),
                  "Mac_moveSteps": np.array([pars["Mac_moveSteps"]])}
    
    print('Param TCell_moveSteps: ')
    print(np.array([pars["TCell_moveSteps"]]))
    print('Param Fib_moveSteps: ')
    print(np.array([pars["Fib_moveSteps"]]))
    print('Param Mac_moveSteps: ')
    print(np.array([pars["Mac_moveSteps"]]))

    # Run the replicates
    #for replicateID in range(model.numReplicates):
    #model.RunModel(SampleID, replicateID, Parameters, RemoveConfigFile=True)
    #output_folders.append(model.get_outputPath(SampleID, replicateID))

    model.RunModel(SampleID, Parameters)
    #output_folders.append(model.get_outputPath())
    output_folders = model.get_outputPath()
    # Return the sum stats of replicates
    return {
           'Reference_CD8T_Weight_Tumor': model_summary(output_folders, 'Reference_CD8T_Weight_Tumor'),
           'Reference_CD4T_Weight_Tumor': model_summary(output_folders, 'Reference_CD4T_Weight_Tumor'),
           'Reference_Fibroblasts_Weight_Fibroblasts': model_summary(output_folders, 'Reference_Fibroblasts_Weight_Fibroblasts'),
           'Reference_CD8T_Weight_Macrophages': model_summary(output_folders, 'Reference_CD8T_Weight_Macrophages'),
           'Reference_Tregs_Weight_Tumor': model_summary(output_folders, 'Reference_Tregs_Weight_Tumor'),
           'Reference_Fibroblasts_Weight_Tregs': model_summary(output_folders, 'Reference_Fibroblasts_Weight_Tregs'),
           'Reference_Tregs_Weight_Tregs': model_summary(output_folders, 'Reference_Tregs_Weight_Tregs')
       }

############ Setup Dask cluster ############
# Each work run a PhysiCell Model

#cluster = SLURMCluster(cores=1, memory="48GB", shebang='#!/bin/bash', queue="parallel", walltime="08:00:00")

#print(cluster.job_script())
#client = Client(cluster)
#client

# Allocate the Total Cores = num_workers * PhysiCellModel.omp_num_threads
#num_workers = 2
#cluster.scale(num_workers) 


############ Calibration ############
# metric function
def euclidean_dist(simulation, data): # distance function
    dist_cd8t = np.linalg.norm( (np.array(data['Reference_CD8T_Weight_Tumor']) - np.array(simulation['Reference_CD8T_Weight_Tumor'])) / np.array(data['Reference_CD8T_Weight_Tumor_sd']))
    dist_cd4t = np.linalg.norm( (np.array(data['Reference_CD4T_Weight_Tumor']) - np.array(simulation['Reference_CD4T_Weight_Tumor'])) / np.array(data['Reference_CD4T_Weight_Tumor_sd']))
    dist_Fib = np.linalg.norm( (np.array(data['Reference_Fibroblasts_Weight_Fibroblasts']) - np.array(simulation['Reference_Fibroblasts_Weight_Fibroblasts'])) / np.array(data['Reference_Fibroblasts_Weight_Fibroblasts_sd']))
    dist_Mac = np.linalg.norm( (np.array(data['Reference_CD8T_Weight_Macrophages']) - np.array(simulation['Reference_CD8T_Weight_Macrophages'])) / np.array(data['Reference_CD8T_Weight_Macrophages_iqr']))
    dist_treg = np.linalg.norm( (np.array(data['Reference_Tregs_Weight_Tumor']) - np.array(simulation['Reference_Tregs_Weight_Tumor'])) / np.array(data['Reference_Tregs_Weight_Tumor_sd']))
    #dist_fib_treg = np.linalg.norm( (np.array(data['Reference_Fibroblasts_Weight_Tregs']) - np.array(simulation['Reference_Fibroblasts_Weight_Tregs'])) / np.array(data['Reference_Fibroblasts_Weight_Tregs_iqr']))
    dist_treg_treg = np.linalg.norm( (np.array(data['Reference_Tregs_Weight_Tregs']) - np.array(simulation['Reference_Tregs_Weight_Tregs'])) / np.array(data['Reference_Tregs_Weight_Tregs_iqr']))
    return dist_cd8t + dist_cd4t + dist_Fib + dist_Mac + dist_treg

# prior distribution
# bounds1 = [0.1, 100.0]; loc1 = 100.0; scale1 = 0.5; 
# bounds2 = [0.01, 1.0]; loc2 = 0.1; scale2 = 0.5; 
parameter_prior = Distribution( TCell_moveSteps=RV("uniform", 20, 40), Fib_moveSteps=RV("uniform", 5, 15), Mac_moveSteps=RV("uniform", 5, 15) )
parameter_prior.get_parameter_names()

# define sampler
fitting_metric = 'all_calibration_8_treg_tumor'
#dask_sampler = sampler.DaskDistributedSampler(client)

#abc = ABCSMC(models=RunSPQSPmodel, parameter_priors=parameter_prior, distance_function=euclidean_dist, population_size=50, sampler=dask_sampler)
abc = ABCSMC(models=RunSPQSPmodel, parameter_priors=parameter_prior, distance_function=euclidean_dist, population_size=100)
# observed data 
theta1_true, theta2_true = np.array([10.0, 1.0]) # drug_dose, dna_damage_rate
measurement_data  = {
       'Reference_CD8T_Weight_Tumor': np.array([0.36]),  
       'Reference_CD4T_Weight_Tumor': np.array([0.43]),  
       'Reference_Fibroblasts_Weight_Fibroblasts': np.array([0.285]), 
       'Reference_CD8T_Weight_Macrophages': np.array([0.07]),  
       'Reference_Tregs_Weight_Tumor': np.array([0.44]),
       'Reference_Fibroblasts_Weight_Tregs': np.array([0.006]),
       'Reference_Tregs_Weight_Tregs': np.array([0.045]),
       'Reference_CD8T_Weight_Tumor_sd': np.array([0.17]),  
       'Reference_CD4T_Weight_Tumor_sd': np.array([0.13]),  
       'Reference_Fibroblasts_Weight_Fibroblasts_sd': np.array([0.14]),
       'Reference_CD8T_Weight_Macrophages_sd': np.array([0.06]),
       'Reference_Tregs_Weight_Tumor_sd': np.array([0.22]),
       'Reference_Fibroblasts_Weight_Tregs_sd': np.array([0.007]),
       'Reference_Tregs_Weight_Tregs_sd': np.array([0.08]),
       'Reference_CD8T_Weight_Tumor_iqr': np.array([0.23]),  
       'Reference_CD4T_Weight_Tumor_iqr': np.array([0.17]),  
       'Reference_Fibroblasts_Weight_Fibroblasts_iqr': np.array([0.21]),  
       'Reference_CD8T_Weight_Macrophages_iqr': np.array([0.08]),
       'Reference_Tregs_Weight_Tumor_iqr': np.array([0.23]),
       'Reference_Fibroblasts_Weight_Tregs_iqr': np.array([0.008]),
       'Reference_Tregs_Weight_Tregs_iqr': np.array([0.06])
   } 

# initialize a new ABC inference run
db_path = "sqlite:///" + os.path.join(f"pyABC_result/{fitting_metric}.db")
abc.new(db_path, measurement_data)

history = abc.run(max_nr_populations=4, minimum_epsilon=0.01)

from pyabc import storage
import matplotlib.pyplot as plt
from pyabc import ABCSMC, Distribution, RV, sampler
from pyabc import LocalTransition, MedianEpsilon
from pyabc.visualization import plot_data_callback, plot_kde_2d, plot_kde_1d
import numpy as np
import math

fitting_metric = 'all_calibration_8_treg_tumor'
history = storage.History(f"sqlite:///pyABC_result/{fitting_metric}.db")

num_plots = history.max_t + 1
num_rows = math.ceil(num_plots / 3)  # Use 3 columns instead of 2
num_cols = min(3, num_plots)  # Maximum of 3 columns

fig = plt.figure(figsize=(10, 12))
for t in range(history.max_t + 1):
    ax = fig.add_subplot(num_rows, num_cols, t + 1)

    ax = plot_kde_1d(
        *history.get_distribution(m=0, t=t),
        'TCell_moveSteps',
        xmin=0,
        xmax=100,
        numx=200,
        ax=ax,
    )
    ax.set_title(f"Posterior t={t}")

    ax.legend()
fig.tight_layout()
plt.savefig(f'fitting_result_1d_TCell_moveSteps.png')


num_plots = history.max_t + 1
num_rows = math.ceil(num_plots / 3)  # Use 3 columns instead of 2
num_cols = min(3, num_plots)  # Maximum of 3 columns


fig2 = plt.figure(figsize=(10, 12))
for t in range(history.max_t + 1):
    ax = fig2.add_subplot(num_rows, num_cols, t + 1)

    ax = plot_kde_2d(
        *history.get_distribution(m=0, t=t),
        "Fib_moveSteps",
        "Mac_moveSteps",
        xmin=0,
        xmax=20,
        numx=200,
        ymin=0,
        ymax=20,
        numy=200,
        ax=ax,
    )
    ax.set_title(f"Posterior t={t}")

    ax.legend()
fig2.tight_layout()
plt.savefig('fitting_result_2d.png')
