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
trigeminal_measures_data_raw = trigeminal_measures_data_raw.rename(columns={
    "R28": "ammo",
    "Lateralisierung (x/20)": "Lateralisierung",
    "CO2-Schwelle": "CO2"
})
ammo = trigeminal_measures_data_raw["ammo"]
ammo_reflect = ammo.max(skipna=True) + 1 - ammo
trigeminal_measures_data_raw["ammo_reflect"] = ammo_reflect
TrigeminalVariableNames = trigeminal_measures_data_raw.columns

for variable in TrigeminalVariableNames:
    data_subset = copy.copy(trigeminal_measures_data_raw[variable])
    if variable == "CO2":
        # 1. All cases (original name)
        explore_tukey_lop(data=data_subset)
        # 2. Cases up to 559, renamed
        data_first559 = data_subset.loc[:548].to_frame().rename(columns={"CO2": "CO2_first549"})
        explore_tukey_lop(data=data_first559["CO2_first549"])
        # 3. Cases from 559 to end, renamed
        data_last559 = data_subset.loc[548:].to_frame().rename(columns={"CO2": "CO2_550toLast"})
        explore_tukey_lop(data=data_last559["CO2_550toLast"])
    else:
        # All cases for other variables
        explore_tukey_lop(data=data_subset)