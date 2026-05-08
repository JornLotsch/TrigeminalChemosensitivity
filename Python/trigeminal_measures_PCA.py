#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 16:14:53 2022

@author: joern
"""

# %% imports
import os
os.chdir("/home/joern/Aktuell/TrigeminalSensitivity/08AnalyseProgramme/Python/")

import pandas as pd

from sklearn.preprocessing import StandardScaler
from compute_PCA_altColors import perform_pca
from cABC_analysis import  cABC_analysis


# %% Functions

# %% PCA

PCA_data = pd.read_csv(
    "/home/joern/Aktuell/TrigeminalSensitivity/08AnalyseProgramme/R/" +  "heat_matrix_all_4_PCA.csv", index_col=0)

PCA_data = pd.DataFrame(StandardScaler().fit_transform(
    PCA_data), columns=PCA_data .columns)

PCA_trigeminal, PCA_trigeminal_feature_importance = perform_pca(
    PCA_data, PC_criterion="KaiserGuttman", plotReduced=3)

PCA_features =  cABC_analysis(PCA_trigeminal_feature_importance)
PCA_features 