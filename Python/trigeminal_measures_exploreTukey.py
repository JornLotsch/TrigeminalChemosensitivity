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
import copy

from explore_tukey_lop import explore_tukey_lop

# %% Read data
pfad_o = "/home/joern/Aktuell/TrigeminalSensitivity/"
pfad_r = "08AnalyseProgramme/R/"
pfad_py = "08AnalyseProgramme/Python/"
pfad_u3 = "09Originale/"

filename = "trigeminal_measures_data.csv"

trigeminal_measures_data_raw = pd.read_csv(pfad_o + pfad_r + filename, index_col=0)
print(trigeminal_measures_data_raw)
AmmoLa_intensity = trigeminal_measures_data_raw["AmmoLa_intensity"]
AmmoLa_intensity_reflect = AmmoLa_intensity.max(skipna=True) + 1 - AmmoLa_intensity
trigeminal_measures_data_raw["AmmoLa_intensity_reflect"] = AmmoLa_intensity_reflect
TrigeminalVariableNames = trigeminal_measures_data_raw.columns

for variable in TrigeminalVariableNames:
    data_subset = copy.copy(trigeminal_measures_data_raw[variable])
    if variable == "CO2_threshold":
        # 1. All cases (original name)
        figure = explore_tukey_lop(data=data_subset)
        # 2. Cases up to 559, renamed
        data_first559 = data_subset.loc[:548].to_frame().rename(columns={"CO2_threshold": "CO2_threshold_first549"})
        figure = explore_tukey_lop(data=data_first559["CO2_threshold_first549"])
        # 3. Cases from 559 to end, renamed
        data_last559 = data_subset.loc[548:].to_frame().rename(columns={"CO2_threshold": "CO2_threshold_550toLast"})
        figure = explore_tukey_lop(data=data_last559["CO2_threshold_550toLast"])
    else:
        # All cases for other variables
        figure = explore_tukey_lop(data=data_subset)
        
# Example code for file saving
# figure = explore_tukey_lop(data_subset, save_fig=True)  # Save to default name
# figure = explore_tukey_lop(data_subset, save_fig=False)  # Do not save
# figure = explore_tukey_lop(data_subset, save_fig=True, fig_path="custom_name.svg")  # Save to custom file

