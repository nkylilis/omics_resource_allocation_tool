#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 16:09:03 2021

@author: nicolaskylilis
"""
# packages
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")



#### analysis results folder
dirname = 'analysis_results'
if os.path.exists(dirname):
   import shutil
   shutil.rmtree(dirname)
   os.mkdir(dirname)
else:
    os.mkdir(dirname)


#%% Mapper: KEGG Orthology to KEGG Pathway 
# Note: https://www.kegg.jp/brite/ko00001
# Hierarchical Classification KEGG BRITE file: 00001 captures the functional relationships of KO ortholog objects (genes and proteins) to Pathway maps objects
# An ortholog maybbe related to more than one Pathway object

#### Load mapper from local
# load KO orthology database mapping to high-level function - species specific
from kegg_tools import kegg_brite
eco_brite       = kegg_brite.load_db_pathways_hierarchy(org_id="eco").set_index("D")
# eco_brite_dupl  = eco_brite[eco_brite.index.duplicated()]


#%% KEGG Genes db

#### Load organsim genes from KEGG genes database (local)
from kegg_tools import kegg_genes
db_genes            = kegg_genes.load_db()
db_genes            = db_genes[db_genes["ENTRY_ORG"] == "T00007"].set_index("ENTRY") # T00007 = e. coli
#duplicates         = db_genes[db_genes.index.duplicated()]


#### Assign KEGG Pathways identifiers to organisms genes
db_genes_pathways   = pd.DataFrame()
for i, row in db_genes.iterrows():
    try:        
        temp    = eco_brite.loc[[i]]
        n       = temp.shape[0]
        for c in range(0,n):
            new_row             = pd.concat([row.to_frame().T, temp.iloc[[c]]], axis="columns")
            db_genes_pathways   = pd.concat([db_genes_pathways, new_row], axis="index")
    except: pass
## filter for protein encoding genes
db_genes_pathways = db_genes_pathways[db_genes_pathways["ENTRY_TYPE"] == "CDS"]


## Keep assignments to KEGG Pathways: accumulated knowledge of cellular and organism-level functions represented in terms of molecular interaction and reaction networks and are categorized into:
# KEGG Pathways
# 9110 Metabolism
# 9120 Genetic Information Processing
# 9130 Environmental Information Processing
# 9140 Cellular Processes
# 9150 Organismal Systems
# 9160 Human Diseases
db_genes_9180       = db_genes_pathways[db_genes_pathways['A'] == 9180]                                         # Brite Hieraerhies
db_genes_9190       = db_genes_pathways[db_genes_pathways['A'] == 9190]                                         # Not Included in Pathway or Brite
db_genes_pathways   = db_genes_pathways[(db_genes_pathways['A'] != 9180) & (db_genes_pathways['A'] != 9190)]    # Pathways [9110-9160]                     


del temp, row, new_row, n, i, c

#%% Proteomics

## load data
df_proteome     = pd.read_csv("proteomics/preprocessed_data/C2.csv", index_col=(0))

## assign mass to genes/pathways ( Note: if a KO belong to multiple pathways, mass is distributed equally)
mass_lst    =[]
for i,row in db_genes_pathways.iterrows():
    try: 
        mass    = df_proteome["mass fraction"].loc[i]
        n       = db_genes_pathways.loc[[i]].shape[0]
        mass_lst+=[mass/n]
    except: mass_lst+=[np.NaN]
db_genes_pathways["mass fraction"] =  mass_lst



## assign mass to Brite Hierarhies
db_genes_brite      = db_genes_9180.loc[db_genes_9180.index.difference(db_genes_pathways.index)]
mass_lst    =[]
for i,row in db_genes_brite.iterrows():
    try: 
        mass    = df_proteome["mass fraction"].loc[i]
        n       = db_genes_brite.loc[[i]].shape[0]
        mass_lst+=[mass/n]
    except: mass_lst+=[np.NaN]
db_genes_brite["mass fraction"] =  mass_lst   


## assign mass to Not Included in Pathway or Brite
mass_lst    =[]
for i,row in db_genes_9190.iterrows():
    try: 
        mass    = df_proteome["mass fraction"].loc[i]
        n       = db_genes_9190.loc[[i]].shape[0]
        mass_lst+=[mass/n]
    except: mass_lst+=[np.NaN]
db_genes_9190["mass fraction"] =  mass_lst


## unmapped to kid
mapped_lst      = list(set(db_genes_pathways.index.tolist() + db_genes_brite.index.tolist() + db_genes_9190.index.tolist()))
unmapped        = df_proteome.loc[df_proteome.index.difference(mapped_lst)]


## save results table
mP = db_genes_pathways["mass fraction"].sum()
mB = db_genes_brite["mass fraction"].sum()
mB = db_genes_brite["mass fraction"].sum()
mN = db_genes_9190["mass fraction"].sum()
mU = unmapped['mass fraction'].sum()

df = pd.DataFrame([mP,mB,mN,mU], columns=['fraction'], index=["KEGG Pathway mapped", 'Brite Hierarhies mapped (not Pathway)', 'kid mapped (not Pathway of Brite)' ,  "unmapped (no kid)"])
df.index.name = 'index'
fpath = os.path.join(dirname,'kid_mapped_proteome.csv')
df.to_csv(fpath, index=True, sep='\t')

    
# Percentage of proteome mapped to KO id
# plot
import matplotlib.pyplot as plt
fig, ax         = plt.subplots(1,1, figsize=(6,4))

x               = [mP, mB, mN, mU]
labels          = ["KEGG Pathway mapped", 'Brite Hierarhies mapped (not Pathway)', 'kid mapped (not Pathway of Brite)' ,  "unmapped (no kid)"]
colors          = ["gold", 'salmon', 'cyan', "grey"]

wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, colors = colors)
ax.set_title("Proteome mapped to KO id")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

fpath = os.path.join(dirname,'kid_mapped_proteome.png')
plt.savefig(fpath)



#%% Resource allocation - Level A

pathways        = sorted(list(set(db_genes_pathways["A"].tolist())))

path_lst =["unmapped (no kid)",'kid mapped (not Pathway of Brite)', 'Brite Hierarhies mapped (not Pathway)']
mass_lst =[mU, mN, mB]
for p in pathways:
    mass    = db_genes_pathways[db_genes_pathways["A"] == p]["mass fraction"].sum()
    path_lst +=[p]
    mass_lst +=[mass]

df_A = pd.DataFrame(zip(path_lst,mass_lst), columns=["Pathway","mass fraction"])

# add pathways names
names_mapper    = kegg_brite.load_db_pathways_names_mapper(org_id='eco')
df_A["name"]    = ["unmapped (no kid)",'kid mapped (not Pathway of Brite)', 'Brite Hierarhies mapped (not Pathway)'] + names_mapper["names"].loc[df_A["Pathway"].tolist()[3:]].tolist()


#### save table to file
fpath = os.path.join(dirname,'Resource_allocation_Level_A.csv')
df_A.to_csv(fpath, index=False, sep='\t')

#### plot
import matplotlib.pyplot as plt
fig, ax         = plt.subplots(1,1, figsize=(12,4))
fig.subplots_adjust(left=-0.5)

x               = df_A["mass fraction"].tolist()
labels          = df_A["name"].tolist()
colors          = ["grey",'cyan','salmon',"gold","orangered","limegreen", "skyblue", "plum", "indigo", "green"]

explode = np.zeros(len(labels))
explode[0:3]=0.2
wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, colors = colors, explode=None)

for w in wedges[0:3]: 
    w.set_alpha(0.2)
    
ax.set_title("Proteome mapped to KEGG Pathways - Level A")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
fpath = os.path.join(dirname,'Resource_allocation_Level_A.png')
plt.savefig(fpath)
    
#%% Resource allocation - Level B


pathways        = sorted(list(set(db_genes_pathways["B"].tolist())))

path_lst =["unmapped (no kid)",'kid mapped (not Pathway of Brite)', 'Brite Hierarhies mapped (not Pathway)']
mass_lst =[mU, mN, mB]
a_lst   = ["unmapped (no kid)",'kid mapped (not Pathway of Brite)', 'Brite Hierarhies mapped (not Pathway)']
for p in pathways:
    temp    = db_genes_pathways[db_genes_pathways["B"] == p]
    mass    = temp["mass fraction"].sum()
    path_lst +=[p]
    mass_lst +=[mass]
    a_lst   += list(set(temp["A"]))
    

df_B = pd.DataFrame(zip(path_lst,mass_lst, a_lst), columns=["Pathway","mass fraction","Level_A"])

# add pathways names
names_mapper    = kegg_brite.load_db_pathways_names_mapper(org_id="eco")
a_lst = ["unmapped (no kid)",'kid mapped (not Pathway of Brite)', 'Brite Hierarhies mapped (not Pathway)'] + names_mapper["names"].loc[df_B["Level_A"].tolist()[3:]].tolist()
b_lst = ['N/A','N/A','N/A'] + names_mapper["names"].loc[df_B["Pathway"].tolist()[3:]].tolist()
names = []
for a,b in zip(a_lst,b_lst):
    name  = a + "  |  " + b
    names += [name]
    
df_B["Pathway name"]    = names


# sort values
df_Btemp        = df_B.iloc[0:3]
df_B            = df_B.iloc[3:]
df_B.sort_values(by=["mass fraction"], ascending=[False], inplace=True)
df_B            = pd.concat([df_Btemp,df_B],axis=0)


#### save table to file
fpath = os.path.join(dirname,'Resource_allocation_Level_B.csv')
df_B.to_csv(fpath, index=False, sep='\t')


#### plot
import matplotlib.pyplot as plt

fig, ax         = plt.subplots(1,1, figsize=(12,4))
fig.subplots_adjust(left=-0.5)

x               = df_B["mass fraction"].iloc[0:12].tolist() + [df_B["mass fraction"].iloc[12:].sum()]
labels          = df_B["Pathway name"].iloc[0:12].tolist() + ['others']
colors          = ["grey",'cyan','salmon',"gold","orangered","limegreen", "skyblue", "plum", "indigo", "green",'yellow','royalblue','lightgrey']
explode = np.zeros(len(labels))
explode[0:3]=0.2
wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, explode=None, colors = colors)

for w in wedges[0:3]: 
    w.set_alpha(0.2)
    
ax.set_title("Proteome mapped to KEGG Pathways - Level B")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
fpath = os.path.join(dirname,'Resource_allocation_Level_B.png')
plt.savefig(fpath)

#%% Resource allocation - Level C

pathways        = sorted(list(set(db_genes_pathways["C"].tolist())))

path_lst =["unmapped (no kid)",'kid mapped (not Pathway of Brite)', 'Brite Hierarhies mapped (not Pathway)']
mass_lst =[mU, mN, mB]
a_lst   = ["unmapped (no kid)",'kid mapped (not Pathway of Brite)', 'Brite Hierarhies mapped (not Pathway)']
b_lst   = ['N/A','N/A','N/A']
c_lst   = ['N/A','N/A','N/A']
for p in pathways:
    temp    = db_genes_pathways[db_genes_pathways["C"] == p]
    mass    = temp["mass fraction"].sum()
    path_lst +=[p]
    mass_lst +=[mass]
    a_lst   += list(set(temp["A"]))
    b_lst   += list(set(temp["B"]))
    c_lst   += list(set(temp["C"]))
    

df_C = pd.DataFrame(zip(path_lst,mass_lst, a_lst, b_lst,c_lst), columns=["Pathway","mass fraction","Level_A",'Level_B','Level_C'])

# add pathways names
names_mapper    = kegg_brite.load_db_pathways_names_mapper(org_id="eco")
a_lst = a_lst[0:3] + names_mapper["names"].loc[df_C["Level_A"].tolist()[3:]].tolist()
b_lst = b_lst[0:3] + names_mapper["names"].loc[df_C["Level_B"].tolist()[3:]].tolist()
c_lst = b_lst[0:3] + names_mapper["names"].loc[df_C["Level_C"].tolist()[3:]].tolist()
names = []
for a,b,c in zip(a_lst,b_lst, c_lst):
    name  = a + "  |  " + b+ "  |  " + c
    names += [name]
    
df_C["Pathway name"]    = names


# sort values
df_Ctemp        = df_C.iloc[0:3]
df_C            = df_C.iloc[3:]
df_C.sort_values(by=["mass fraction"], ascending=[False], inplace=True)
df_C            = pd.concat([df_Ctemp,df_C],axis=0)



#### save table to file
fpath = os.path.join(dirname,'Resource_allocation_Level_C.csv')
df_C.to_csv(fpath, index=False, sep='\t')


#### plot
import matplotlib.pyplot as plt

fig, ax         = plt.subplots(1,1, figsize=(12,4))
fig.subplots_adjust(left=-0.5)


x               = df_C["mass fraction"].iloc[0:12].tolist() + [df_C["mass fraction"].iloc[12:].sum()]
labels          = df_C["Pathway name"].iloc[0:12].tolist() + ['others']
colors          = ["grey",'cyan','salmon',"gold","orangered","limegreen", "skyblue", "plum", "indigo", "green",'yellow','royalblue','grey']

# explode = np.zeros(len(labels))
# explode[3:]=0.5

wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, colors = colors) #explode=explode,
for w in wedges[0:3]: 
    w.set_alpha(0.2)
    
ax.set_title("Proteome mapped to KEGG Pathways - Level B")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
fpath = os.path.join(dirname,'Resource_allocation_Level_C.png')
plt.savefig(fpath)
