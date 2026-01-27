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

filename = "training_data_to_transform.csv"

trigeminal_measures_data_raw = pd.read_csv(pfad_o + pfad_r + filename, index_col=0)
print(trigeminal_measures_data_raw)
AmmoLa_intensity = trigeminal_measures_data_raw["AmmoLa intensity"]
AmmoLa_intensity_reflect = AmmoLa_intensity.max(skipna=True) + 1 - AmmoLa_intensity
CO2_threshold = trigeminal_measures_data_raw["CO2 threshold"]
CO2_threshold_reflect = CO2_threshold.max(skipna=True) + 1 - CO2_threshold
Lateralization = trigeminal_measures_data_raw["Lateralization (x/20)"]
Lateralization_reflect = Lateralization.max(skipna=True) + 1 - Lateralization
trigeminal_measures_data_raw["AmmoLa_intensity_reflect"] = AmmoLa_intensity_reflect
trigeminal_measures_data_raw["CO2_threshold_reflect"] = CO2_threshold_reflect
trigeminal_measures_data_raw["Lateralization_reflect"] = Lateralization_reflect


TrigeminalVariableNames = trigeminal_measures_data_raw.columns

used_powers = [-2, -1, -0.5, 0, 0.5, 1, 2]

for variable in TrigeminalVariableNames:
    data_subset = copy.copy(trigeminal_measures_data_raw[variable])
    if variable == "CO2 threshold":
        # 1. All cases (original name)
        figure = explore_tukey_lop(data=data_subset, powers=used_powers, max_points_plot=500)
        # 2. Cases up to 559, renamed
        data_first559 = data_subset.loc[:548].to_frame().rename(columns={"CO2 threshold": "CO2_threshold_first549"})
        figure = explore_tukey_lop(data=data_first559["CO2_threshold_first549"], powers=used_powers, max_points_plot=500)
        # 3. Cases from 559 to end, renamed
        data_last559 = data_subset.loc[548:].to_frame().rename(columns={"CO2 threshold": "CO2_threshold_550toLast"})
        figure = explore_tukey_lop(data=data_last559["CO2_threshold_550toLast"], powers=used_powers, max_points_plot=500)
    elif variable == "CO2_threshold_reflect":
        # 1. All cases (original name)
        figure = explore_tukey_lop(data=data_subset, powers=used_powers, max_points_plot=500)
        # 2. Cases up to 559, renamed
        data_first559 = data_subset.loc[:548].to_frame().rename(columns={"CO2_threshold_reflect": "CO2_threshold_reflect_first549"})
        figure = explore_tukey_lop(data=data_first559["CO2_threshold_reflect_first549"], powers=used_powers, max_points_plot=500)
        # 3. Cases from 559 to end, renamed
        data_last559 = data_subset.loc[548:].to_frame().rename(columns={"CO2_threshold_reflect": "CO2_threshold_reflect_550toLast"})
        figure = explore_tukey_lop(data=data_last559["CO2_threshold_reflect_550toLast"], powers=used_powers, max_points_plot=500)
    else:
        # FIXED: Pass powers=used_powers for ALL variables
        figure = explore_tukey_lop(data=data_subset, powers=used_powers, max_points_plot=500)

        
# File saving
# figure = explore_tukey_lop(data_subset, save_fig=True)  # Save to default name
# figure = explore_tukey_lop(data_subset, save_fig=False)  # Do not save
# figure = explore_tukey_lop(data_subset, save_fig=True, fig_path="custom_name.svg")  # Save to custom file

