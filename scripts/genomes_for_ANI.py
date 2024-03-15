#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 11:05:52 2024

Copy all phasco genomes into the folder

@author: ekateria
"""

import pandas as pd
import shutil
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import numpy as np

genomes_path='genomes_names_phasco.tsv'

genomes=pd.read_csv(genomes_path, sep='\t')
genomes.columns=['Path','Species']
genomes['Sample']=genomes['Path'].apply(lambda row: row.split('_')[0])
genomes['NewName']=genomes['Species']+'_'+genomes['Sample']+'.fasta'

source_path='work_dir/'
dest_path='Genomes_for_ANI/'

for _,file in genomes.iterrows():
    try:
        shutil.copy(source_path+file.Sample+'/binning/DASTool/bins/'+file.Path+'.fasta', dest_path+file.NewName)
    except IOError as e:
        print("Unable to copy file. %s" % e)
        

##This part of the script is performed after the ANI was calculated outside of the script using 
#average_nucleotide_identity.py -o /PyANI_PhascoCRCbiome_ANI_tetra/
# -i /Genomes_for_ANI/ -m TETRA
       
ani=pd.read_csv('ANI_correlations.tab', sep='\t')
ani=ani.set_index('Unnamed: 0')

ani_dissim=1-ani
mds=MDS(n_components=2,dissimilarity='precomputed',random_state=399)
embedding = mds.fit_transform(ani_dissim)


fig=plt.scatter(embedding[0:130, 0], embedding[0:130, 1], c='#3c5488ff') #Pfaecium
plt.scatter(embedding[131:163, 0], embedding[131:163, 1], c='#f39b7f') #Pspp
plt.scatter(embedding[164:, 0], embedding[164:, 1], c='#00a087ff') #Psuccin
plt.plot([-0.24,0.11],[0,0],linestyle='--',color='gray',linewidth=1)
plt.plot([0,0],[-0.17,0.07],linestyle='--',color='gray',linewidth=1)
plt.ylim([-0.17,0.07])
plt.xlim([-0.24,0.11])

#Calculate average within-species ANI

species=['Pfaecium','Pspecies','Psuccinatutens']

def ani_within(sp,data):
    spdata=data.filter(like=sp,axis=0).filter(like=sp,axis=1)
    avg=spdata.values.mean()
    std=spdata.values.std()
    
    return avg, std

for sp in species:
    avg, std=ani_within(sp,ani)
    print(sp+' average ANI: ' + str(avg)+'±' + str(std))

def ani_between(sp,data):
    

    avg=spdata.values.mean()
    std=spdata.values.std()

    return avg, std

Vecs=[]
for sp in species:
     
    spdata=ani.filter(like=sp,axis=0)
    spdata=spdata.loc[:,~spdata.columns.str.contains(sp)]
    
    vec=spdata.values.flatten()
    Vecs.extend(vec)

avg_between=np.mean(Vecs)
std_between=np.std(Vecs)

print('Between species average ANI: ' + str(avg_between)+'±' + str(std_between))
