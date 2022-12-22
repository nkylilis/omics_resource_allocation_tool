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
db_genes_brite   = pd.DataFrame()
for i, row in db_genes.iterrows():
    try:        
        temp    = eco_brite.loc[[i]]
        n       = temp.shape[0]
        for c in range(0,n):
            new_row             = pd.concat([row.to_frame().T, temp.iloc[[c]]], axis="columns")
            db_genes_brite      = pd.concat([db_genes_brite, new_row], axis="index")
    except: pass
## filter for protein encoding genes
# db_genes_brite = db_genes_brite[db_genes_brite["ENTRY_TYPE"] == "CDS"]


## KEGG Pathways: accumulated knowledge of cellular and organism-level functions represented in terms of molecular interaction and reaction networks and are categorized into:
# 9100 Metabolism
# 9120 Genetic Information Processing
# 9130 Environmental Information Processing
# 9140 Cellular Processes
# 9150 Organismal Systems
# 9160 Human Diseases                   

db_genes_pathways_ext   = db_genes_brite.loc[db_genes_brite['A'].isin([9100,9120,9130,9140,9150,9160])] 
db_genes_pathways       = db_genes_brite.loc[db_genes_brite['A'].isin([9100,9120,9130,9140])]
db_genes_pathways_not   = db_genes_brite.loc[db_genes_brite.index.difference(db_genes_pathways_ext.index)]

# Brite Hieraerhies
# manual curation of ontology:  
db_genes_brite_m        = db_genes_brite[db_genes_brite['C'] != 3029] # 03029 Mitochondrial Biogenesis (no mitochondria in E. coli)
# test = db_genes_brite.loc[db_genes_brite.index.difference(db_genes_brite_m.index)] # maing sure I do not lose any genes
db_genes_9180           = db_genes_brite_m[db_genes_brite_m['A'] == 9180]
db_genes_9180_not       = db_genes_brite_m.loc[db_genes_brite_m.index.difference(db_genes_9180.index)]


del temp, row, new_row, n, i, c


#%% Proteomics

## load data
df_proteome     = pd.read_csv("proteomics/preprocessed_data/C2.csv", index_col=(0))


#### KEGG Pathways 
# assign mass to genes/pathways ( Note: if a KO belong to multiple pathways, mass is distributed equally)
mass_lst    =[]
for i,row in db_genes_pathways.iterrows():
    try: 
        mass    = df_proteome["mass fraction"].loc[i]
        n       = db_genes_pathways.loc[[i]].shape[0]
        mass_lst+=[mass/n]
    except: mass_lst+=[np.NaN]
db_genes_pathways["mass fraction"] =  mass_lst

## KEGG Pathways extended
# assign mass to genes/pathways ( Note: if a KO belong to multiple pathways, mass is distributed equally)
mass_lst    =[]
for i,row in db_genes_pathways_ext.iterrows():
    try: 
        mass    = df_proteome["mass fraction"].loc[i]
        n       = db_genes_pathways_ext.loc[[i]].shape[0]
        mass_lst+=[mass/n]
    except: mass_lst+=[np.NaN]
db_genes_pathways_ext["mass fraction"] =  mass_lst

## NOT in KEGG Pathways 
# assign mass to genes/pathways ( Note: if a KO belong to multiple pathways, mass is distributed equally)
mass_lst    =[]
for i,row in db_genes_pathways_not.iterrows():
    try: 
        mass    = df_proteome["mass fraction"].loc[i]
        n       = db_genes_pathways_not.loc[[i]].shape[0]
        mass_lst+=[mass/n]
    except: mass_lst+=[np.NaN]
db_genes_pathways_not["mass fraction"] =  mass_lst


## unmapped to kid
mapped_lst      = list(set(db_genes_pathways_ext.index.tolist() + db_genes_pathways_not.index.tolist()))
unmapped        = df_proteome.loc[df_proteome.index.difference(mapped_lst)]



mP      = db_genes_pathways["mass fraction"].sum()
mPx     = db_genes_pathways_ext["mass fraction"].sum()
mPx_not = db_genes_pathways_not["mass fraction"].sum()
mU      = unmapped['mass fraction'].sum()


# plot
import matplotlib.pyplot as plt
fig, ax         = plt.subplots(1,1, figsize=(12,4))

x               = [mPx, mPx_not, mU]
labels          = ["KEGG Pathway mapped", 'kid mapped (not Pathway)',  "kid unmapped"]
colors          = ["palegreen", 'lightcoral', "grey"]

wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, colors = colors)
ax.set_title("Proteome mapped to KEGG Pathways")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

fpath = os.path.join(dirname,'Pathways_mapped_proteome.png')
plt.savefig(fpath)





#### Brite Hierarhies  
## assign mass to Brite Hierarhies
mass_lst    =[]
for i,row in db_genes_9180.iterrows():
    try: 
        mass    = df_proteome["mass fraction"].loc[i]
        n       = db_genes_9180.loc[[i]].shape[0]
        mass_lst+=[mass/n]
    except: mass_lst+=[np.NaN]
db_genes_9180["mass fraction"] =  mass_lst   


## assign mass to Not  Brite
mass_lst    =[]
for i,row in db_genes_9180_not.iterrows():
    try: 
        mass    = df_proteome["mass fraction"].loc[i]
        n       = db_genes_9180_not.loc[[i]].shape[0]
        mass_lst+=[mass/n]
    except: mass_lst+=[np.NaN]
db_genes_9180_not["mass fraction"] =  mass_lst   


## unmapped to kid
mapped_lst      = list(set(db_genes_9180.index.tolist() + db_genes_9180_not.index.tolist()))
unmapped_b        = df_proteome.loc[df_proteome.index.difference(mapped_lst)]


mB      = db_genes_9180["mass fraction"].sum()
mB_not  = db_genes_9180_not["mass fraction"].sum()
mUb     = unmapped_b['mass fraction'].sum()


# plot
import matplotlib.pyplot as plt
fig, ax         = plt.subplots(1,1, figsize=(12,4))

x               = [mB, mB_not, mUb]
labels          = ['Brite Hierarhies mapped', 'kid mapped (not  Brite)' ,  "kid unmapped "]
colors          = ["palegreen", 'lightcoral', 'grey']

wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, colors = colors)
ax.set_title("Proteome mapped to Brite Hierarhies")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))

fpath = os.path.join(dirname,'Brite_Hierrarhies_mapped_proteome.png')
plt.savefig(fpath)




## save results table
df = pd.DataFrame([mP, mPx, mPx_not, mU, mB, mB_not, mUb], columns=['fraction'], index=["KEGG Pathways", "KEGG Pathways extended", "kid NOT KEGG Pathways", 'Unmapped' , 'Brite Hierarhies', 'kid NOT Brite Hierarhies' ,  "Unmapped b"])
df.index.name = 'index'
fpath = os.path.join(dirname,'kid_mapped_proteome.csv')
df.to_csv(fpath, index=True, sep='\t')

#%% Resource allocation - Level A

color_sch       = dict(zip(sorted(list(set(db_genes_pathways['A'].tolist()))),['yellow','red','limegreen','plum','skyblue']))
pathways        = sorted(list(set(db_genes_pathways["A"].tolist())))

path_lst =["unmapped (no kid)",'kid mapped (not Pathway)']
mass_lst =[mU, mPx_not]
color   = ['black', 'grey']
for p in pathways:
    temp    = db_genes_pathways[db_genes_pathways["A"] == p]
    mass    = temp["mass fraction"].sum()
    path_lst +=[p]
    mass_lst +=[mass]
    color   += [color_sch[list(set(temp["A"]))[0]]]
    
df_A = pd.DataFrame(zip(path_lst,mass_lst), columns=["Pathway","mass fraction"])

# add pathways names
names_mapper    = kegg_brite.load_db_pathways_names_mapper(org_id='eco')
df_A["name"]    = ["unmapped (no kid)",'kid mapped (not Pathway)'] + names_mapper["names"].loc[df_A["Pathway"].tolist()[2:]].tolist()
df_A["colors"]          = color

#### save table to file
fpath = os.path.join(dirname,'Resource_allocation_Level_A.csv')
df_A.to_csv(fpath, index=False, sep='\t')

#### plot
import matplotlib.pyplot as plt
fig, ax         = plt.subplots(1,1, figsize=(12,4))
fig.subplots_adjust(left=-0.5)

x               = df_A["mass fraction"].tolist()
labels          = df_A["name"].tolist()
colors          = df_A['colors'].iloc[0:12].tolist() + ['blueviolet']

explode = np.zeros(len(labels))
explode[0:2]=0.2
wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, colors = colors, explode=None, pctdistance=1.2, wedgeprops= {"edgecolor":"black",'linewidth': 0.25,'antialiased': True})

for w in wedges[0:2]: 
    w.set_alpha(0.2)
    
ax.set_title("Proteome mapped to KEGG Pathways - Level A")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
fpath = os.path.join(dirname,'Resource_allocation_Level_A.png')
plt.savefig(fpath)
    
#%% Resource allocation - KEGG Pathways Level B

color_sch       = dict(zip(sorted(list(set(db_genes_pathways['A'].tolist()))),['yellow','red','limegreen','plum','skyblue']))
pathways        = sorted(list(set(db_genes_pathways["B"].tolist())))

path_lst =["unmapped (no kid)",'kid mapped (not Pathway)']
mass_lst =[mU, mPx_not]
a_lst   = ["unmapped (no kid)",'kid mapped (not Pathway)']
color   = ['black', 'grey']
for p in pathways:
    temp    = db_genes_pathways[db_genes_pathways["B"] == p]
    mass    = temp["mass fraction"].sum()
    path_lst +=[p]
    mass_lst +=[mass]
    a_lst   += list(set(temp["A"]))
    color   += [color_sch[list(set(temp["A"]))[0]]]
    

df_B = pd.DataFrame(zip(path_lst,mass_lst, a_lst), columns=["Pathway","mass fraction","Level_A"])

# add pathways names
names_mapper    = kegg_brite.load_db_pathways_names_mapper(org_id="eco")
a_lst = ["unmapped (no kid)",'kid mapped (not Pathway)'] + names_mapper["names"].loc[df_B["Level_A"].tolist()[2:]].tolist()
b_lst = ['N/A','N/A'] + names_mapper["names"].loc[df_B["Pathway"].tolist()[2:]].tolist()
names = []
for a,b in zip(a_lst,b_lst):
    name  = a + "  |  " + b
    names += [name]
    
df_B["Pathway name"]    = names
df_B["colors"]          = color

# sort values
df_Btemp        = df_B.iloc[0:2]
df_B            = df_B.iloc[2:]
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
colors          = df_B['colors'].iloc[0:12].tolist() + ['blueviolet']
explode = np.zeros(len(labels))
explode[0:2]=0.2
wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, explode=None, colors = colors, pctdistance=1.2,  wedgeprops= {"edgecolor":"black",'linewidth': 0.25,'antialiased': True})

for w in wedges[0:2]: 
    w.set_alpha(0.2)
    
ax.set_title("Proteome mapped to KEGG Pathways - Level B")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
fpath = os.path.join(dirname,'Resource_allocation_Level_B.png')
plt.savefig(fpath)




#%% Resource allocation - KEGG Pathways Level C
color_sch       = dict(zip(sorted(list(set(db_genes_pathways['A'].tolist()))),['yellow','red','limegreen','plum','skyblue']))
pathways        = sorted(list(set(db_genes_pathways["C"].tolist())))

path_lst =["unmapped (no kid)",'kid mapped (not Pathway)']
mass_lst =[mU, mPx_not]
a_lst   = ["unmapped (no kid)",'kid mapped (not Pathway)']
b_lst   = ['N/A','N/A']
c_lst   = ['N/A','N/A']
color   = ['black', 'grey']
for p in pathways:
    temp    = db_genes_pathways[db_genes_pathways["C"] == p]
    mass    = temp["mass fraction"].sum()
    path_lst +=[p]
    mass_lst +=[mass]
    a_lst   += list(set(temp["A"]))
    b_lst   += list(set(temp["B"]))
    c_lst   += list(set(temp["C"]))
    color   += [color_sch[list(set(temp["A"]))[0]]]
    

df_C = pd.DataFrame(zip(path_lst,mass_lst, a_lst, b_lst,c_lst), columns=["Pathway","mass fraction","Level_A",'Level_B','Level_C'])

# add pathways names
names_mapper    = kegg_brite.load_db_pathways_names_mapper(org_id="eco")
a_lst = a_lst[0:2] + names_mapper["names"].loc[df_C["Level_A"].tolist()[2:]].tolist()
b_lst = b_lst[0:2] + names_mapper["names"].loc[df_C["Level_B"].tolist()[2:]].tolist()
c_lst = b_lst[0:2] + names_mapper["names"].loc[df_C["Level_C"].tolist()[2:]].tolist()
names = []
for a,b,c in zip(a_lst,b_lst, c_lst):
    name  = a + "  |  " + b+ "  |  " + c
    names += [name]
    
df_C["Pathway name"]    = names
df_C["colors"]          = color


# sort values
df_Ctemp        = df_C.iloc[0:2]
df_C            = df_C.iloc[2:]
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
colors          = df_C['colors'].iloc[0:12].tolist() + ['blueviolet']
# explode = np.zeros(len(labels))
# explode[3:]=0.5

wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, colors = colors, explode=None, pctdistance=1.2,  wedgeprops= {"edgecolor":"black",'linewidth': 0.25,'antialiased': True})
for w in wedges[0:2]: 
    w.set_alpha(0.2)
    
ax.set_title("Proteome mapped to KEGG Pathways - Level C")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
fpath = os.path.join(dirname,'Resource_allocation_Level_C.png')
plt.savefig(fpath)



#%% Resource allocation - Brite Hierarhies Level B

color_sch    = dict(zip(sorted(list(set(db_genes_9180["B"].tolist()))),['yellow','red','limegreen','plum','skyblue']))
Brite        = sorted(list(set(db_genes_9180["B"].tolist())))

path_lst =["unmapped (no kid)",'kid mapped (not Brite)']
mass_lst =[mUb, mB_not]
color      =["black",'grey']
for p in Brite:
    temp    = db_genes_9180[db_genes_9180["B"] == p]
    mass    = temp["mass fraction"].sum()
    path_lst +=[p]
    mass_lst +=[mass]
    b_lst   += list(set(temp["B"]))
    color   += [color_sch[list(set(temp["B"]))[0]]]

    

df_Br = pd.DataFrame(zip(path_lst,mass_lst, b_lst), columns=["Category","mass fraction","Level B" ])

# add pathways names
names_mapper    = kegg_brite.load_db_pathways_names_mapper(org_id="eco")
b_lst = ["unmapped (no kid)",'kid mapped (not Brite)'] + names_mapper["names"].loc[df_Br["Category"].tolist()[2:]].tolist()
names = []
for b in b_lst:
    name  =  b
    names += [name]
    
df_Br["Pathway name"]    = names
df_Br["colors"]          = color

# sort values
df_Brtemp        = df_Br.iloc[0:2]
df_Br            = df_Br.iloc[2:]
df_Br.sort_values(by=["mass fraction"], ascending=[False], inplace=True)
df_Br            = pd.concat([df_Brtemp,df_Br],axis=0)


#### save table to file
fpath = os.path.join(dirname,'Resource_allocation_BriteHierarhies_Level_B.csv')
df_Br.to_csv(fpath, index=False, sep='\t')


#### plot
import matplotlib.pyplot as plt

fig, ax         = plt.subplots(1,1, figsize=(12,4))
fig.subplots_adjust(left=-0.5)

x               = df_Br["mass fraction"].iloc[0:12].tolist() + [df_Br["mass fraction"].iloc[12:].sum()]
labels          = df_Br["Pathway name"].iloc[0:12].tolist() + ['others']
colors          = df_Br['colors'].iloc[0:12].tolist() + ['blueviolet']
explode = np.zeros(len(labels))
explode[0:2]=0.2
wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, explode=None, colors = colors, pctdistance=1.2,  wedgeprops= {"edgecolor":"black",'linewidth': 0.25,'antialiased': True})

for w in wedges[0:2]: 
    w.set_alpha(0.2)
    
ax.set_title("Proteome mapped to  Brite Hierarhies - Level B")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
fpath = os.path.join(dirname,'Resource_allocation_BriteHierarhies_Level_B.png')
plt.savefig(fpath)


#%% Resource allocation - Brite Hierarhies Level C

color_sch    = dict(zip(sorted(list(set(db_genes_9180["B"].tolist()))),['yellow','red','limegreen','plum','skyblue']))
Brite        = sorted(list(set(db_genes_9180["C"].tolist())))

path_lst    =["unmapped (no kid)",'kid mapped (not Brite)']
mass_lst    =[mUb, mB_not]
b_lst    =["unmapped (no kid)",'kid mapped (not Brite)']
c_lst    =["unmapped (no kid)",'kid mapped (not Brite)']
color      =["black",'grey']
for p in Brite:
    temp    = db_genes_9180[db_genes_9180["C"] == p]
    mass    = temp["mass fraction"].sum()
    path_lst +=[p]
    mass_lst +=[mass]
    b_lst   += list(set(temp["B"]))
    c_lst   += list(set(temp["C"]))
    color   += [color_sch[list(set(temp["B"]))[0]]]

df_Br = pd.DataFrame(zip(path_lst,mass_lst, b_lst, c_lst), columns=["Category","mass fraction","Level B","Level C" ])


# add pathways names
names_mapper    = kegg_brite.load_db_pathways_names_mapper(org_id="eco")
b_lst = ["unmapped (no kid)",'kid mapped (not Brite)'] + names_mapper["names"].loc[df_Br["Level B"].tolist()[2:]].tolist()
c_lst = ['N/A','N/A'] + names_mapper["names"].loc[df_Br["Level C"].tolist()[2:]].tolist()
names = []
for b,c in zip(b_lst, c_lst):
    name  = b+ "  |  " + c
    names += [name]
    
df_Br["Pathway name"]    = names
df_Br["colors"]          = color


# sort values
df_Brtemp        = df_Br.iloc[0:2]
df_Br            = df_Br.iloc[2:]
df_Br.sort_values(by=["mass fraction"], ascending=[False], inplace=True)
df_Br            = pd.concat([df_Brtemp,df_Br],axis=0)


#### save table to file
fpath = os.path.join(dirname,'Resource_allocation_BriteHierarhies_Level_C.csv')
df_Br.to_csv(fpath, index=False, sep='\t')


#### plot
import matplotlib.pyplot as plt

fig, ax         = plt.subplots(1,1, figsize=(12,4))
fig.subplots_adjust(left=-0.5)

x               = df_Br["mass fraction"].iloc[0:12].tolist() + [df_Br["mass fraction"].iloc[12:].sum()]
labels          = df_Br["Pathway name"].iloc[0:12].tolist() + ['others']
colors          = df_Br['colors'].iloc[0:12].tolist() + ['blueviolet']
explode[0:2]=0.2
wedges, texts, autotexts = ax.pie(x, autopct='%1.1f%%', shadow=False, explode=None, colors = colors, pctdistance=1.2,  wedgeprops= {"edgecolor":"black",'linewidth': 0.25,'antialiased': True})

for w in wedges[0:2]: 
    w.set_alpha(0.2)
    
ax.set_title("Proteome mapped to Brite Hierarhies - Level C")
ax.legend(labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
fpath = os.path.join(dirname,'Resource_allocation_BriteHierarhies_Level_C.png')
plt.savefig(fpath)