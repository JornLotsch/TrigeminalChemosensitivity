#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 09:12:28 2023

@author: joern
"""

# %% imports

import os

os.chdir("/home/joern/Aktuell/TrigeminalSensitivity/08AnalyseProgramme/Python/")

import pandas as pd
import numpy as np
import copy
from sklearn.preprocessing import StandardScaler

from explore_tukey_lop import explore_tukey_lop

# %% Read data
pfad_o = "/home/joern/Aktuell/TrigeminalSensitivity/"
pfad_r = "08AnalyseProgramme/R/"
pfad_py = "08AnalyseProgramme/Python/"
pfad_u3 = "09Originale/"

filename = "trigeminal_measures_data.csv"

trigeminal_measures_data_raw = pd.read_csv(pfad_o + pfad_r + filename, index_col=0)
print(trigeminal_measures_data_raw )

trigeminal_measures_data_raw = trigeminal_measures_data_raw.rename(columns={
    "R28": "ammo",
    "Lateralisierung (x/20)": "Lateralisierung",
    "CO2-Schwelle": "CO2"
})


ammo = trigeminal_measures_data_raw["ammo"]
ammo_reflect = ammo.max(skipna=True) + 1 - ammo

trigeminal_measures_data_raw["ammo_reflect"] = ammo_reflect


TrigeminalVariableNames = trigeminal_measures_data_raw.columns

for i, variable in enumerate(TrigeminalVariableNames):
    data_subset = copy.copy(trigeminal_measures_data_raw[variable])
    explore_tukey_lop(data=data_subset)
    
    
