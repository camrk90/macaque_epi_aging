#Import libraries
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import LeaveOneOut
import os

#Python indexes from 0 so call is SAMP - 1
#SAMP=int(os.getenv("SLURM_ARRAY_TASK_ID")) - 1

#Set path
location=input('Are you working locally or remotely today?')

if location=='locally':
    path="/Users/cameronkelsey/Documents/smack_lab/cayo_data/"
elif location=='remotely':
    path="/scratch/ckelsey4/Cayo_meth/"

#Import data
meta=pd.read_table(path + "blood_metadata_full.txt", sep=r'\s+')
pmeth=pd.read_table(path + "pmeth_clock.txt", sep=r'\s+')





