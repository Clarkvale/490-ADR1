#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 14:16:46 2020

@author: benjamin
"""

#ortho database analysis

from databases.orthdb.orthDB import OrthDB
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd



def num_of_orthogroups():
    db = OrthDB()
    sp_count = {}
    species = {'Spar': 62, 'Ctro': 30, 'Smik': 66, 'Spom': 7, 'Cgui': 20, 'Soct': 10, 
               'Clus': 35, 'Sklu': 29, 'Ncra': 57, 'Ylip': 34, 'Calb': 45, 
               'Lelo': 25, 'Agos': 66, 'Klac': 35, 
               'Sbay': 52, 'Anid': 44, 'Dhan': 16, 'Scas': 43, 'Sjap': 30, 
               'Cpar': 18, 'Scer': 199, 'Cgla': 37, 'Kwal': 30}
    
    for s in species.keys():
        db.cursor.execute("SELECT Orthogroup From Promoters WHERE Species = ?", (s,))
        sp_count[s] = len(db.cursor.fetchall())
    db.close()
    return sp_count

def makebar(species_count):
    plt.xticks(rotation = 70)
    bar = sns.barplot(x = list(species_count.keys()), 
                y = list(species_count.values(),))
    bar.figure.savefig("results/num_of_valid_orthogroups_per_species.png")
    