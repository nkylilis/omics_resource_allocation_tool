#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 11:39:59 2022

@author: nicolaskylilis
"""


#### Load datasets
import pandas as pd

df_proteome = pd.read_excel("msb20209536-sup-0010-datasetev9.xlsx")
columns = ['Gene name', 'Gene locus', 'Protein ID', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'D6','D7', 'D8', 'F4', 'F5', 'F6', 'F7', 'F8',"D1",	"D2", "D3",	"D4", "D5", "F2", "F3"]
df_proteome = df_proteome[columns].copy()
df_proteome.set_index("Gene locus", inplace=True)
# look for duplicated gene locus
duplicated = df_proteome[df_proteome.index.duplicated()]
df_proteome = df_proteome[~df_proteome.index.duplicated(keep='first')]


#### Proteome Sectors allocations
# iterate media conditions
conditions = columns[3:]
df_results_proteome = pd.DataFrame([])

for cond in conditions:
    
    # map proteome to protein info
    df_proteome_cond =   df_proteome[[ cond]].copy()
    df_proteome_cond.rename(columns={cond:"mass fraction"}, inplace=True)
    df_proteome_cond.to_csv("preprocessed_data/" + cond + ".csv")
    
