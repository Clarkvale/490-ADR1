#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:58:31 2020

@author: benjamin
"""

from databases.orthdb.orthDB import OrthDB
from fimo_parser.GFF import GFF_Parse
from interval_tree import Motif_Node
from dyad_search import UAS
from Bio.Seq import Seq
import time
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
#from plotly.graph_objects import *
#import plotly.figure_factory as ff 
#from plotly.graph_objects.layout import Layout

def main():
    ##get gff fimo motif hits
    fimo_file = "wg_fimo/ortho_fimo_0.01p.gff"
    print("extracting fimo hits....")
    ts = time.time()
    fimo = GFF_Parse(fimo_file)
    te= time.time()
    print("done in " + str(te-ts) + "s")
    df_fimo = fimo.dataframe
    
    ##pull sequences by index, if they are less than 47bp just take the whole 
    #thing
    ts = time.time()
    print("getting sequences.....")
    seq40 = []
    
    #initialize db connection
    db = OrthDB()
    #pulling sequences by index, if you change the order of the database you have 
    #to re-run everything
    total_scanned_sequence_length = 0
    seq40 = []
    for name in list(set(df_fimo["name"])):
        index = int(name.split("_")[0]) - 1
        
        seq = db.lookup_index(index)[2]
        nodes = []
        for i, row in df_fimo[df_fimo["name"] == name].iterrows():
            if len(seq) <= 47:
                nodes.append(Motif_Node(0,len(seq)))
                break
            elif row["start"] < 20:
                nodes.append(Motif_Node(0, row["end"] + 20))
            elif row["end"] >= len(seq) - 20:
                nodes.append(Motif_Node(row["start"], len(seq)-1))
            else:
                nodes.append(Motif_Node(row["start"] - 20, row["end"] + 20))
        total_promoter_seq = Motif_Node.decompose(nodes)
        total_scanned_sequence_length += sum([node.e - node.s for node in total_promoter_seq])
        for coors in total_promoter_seq:
            seq40.append((name,seq[coors.s:coors.e]))
    te = time.time()
    print("done in " + str(te-ts) + "s")
    
    #look for palindromes from each sequence
    genes = _find_pals(seq40)
    #export everything
    df = export(genes)
    #build heatmap
    return df
def _find_pals(seq40):
    print("finding palindromes......")
    ts = time.time()
    #To find the frequency of each residue
    averages = {"C":[], "G":[], "T":[], "A":[]}
    unique_palindrome_count = 0
    #genes holds all genes with a promoter containing at least one palindrome with
    # a decently high cg content to avoid picking up TATA boxes
    genes = []
    for name, seq in seq40:
        #handling odd sequence length
        if len(seq40) % 2 != 0:
            uas = UAS(name, seq[:-1], 2)
        else:
            uas = UAS(name, seq, 2)
        for key in averages.keys():    
            averages[key].append(seq.count(key)/len(seq))
        if any(uas.pals_in_line):
            try:
                #Overlapping palindromes are collapsed into one set of coordinates
                #and then counted 
                noded_palindromes = [Motif_Node(pal.pos, pal.pos + len(pal)) 
                for pal in uas.pals_in_line]
                unique_palindrome_count += len(Motif_Node.decompose(noded_palindromes))
                
                if(uas.get_best_pal().cg_content() >= 35):
                    genes.append((uas, uas.get_best_pal()))
            except AttributeError:
                print(uas)
                pass
    te = time.time()
    print("done in " + str(te-ts) + "s")
    return genes

def export(genes):
    n = set([x.name for x,y in genes])
    pal_dict = {}
    csv_out = []
    for name,p in genes:
        pal_dict["Species"] = name.name.split("_")[1].split("|")[0]
        pal_dict["Gene"] = name.name.split("|")[1]    
        pal_dict["Sequence"] = p.pal
        pal_dict["Length"] = len(p.pal)
        pal_dict["Score"] = p.score
        csv_out.append(pal_dict.copy())
        pal_dict.clear()
        
    csv = pd.DataFrame(csv_out).drop_duplicates()
    
    csv.to_csv("ortho_palindromes_out.csv")
    
    with open("ortho_uas_gene_names.txt", "w") as fo:
        for name in n:
            fo.write(name + "\n")
    return csv
            
    
    #build export files and heatmap. CSV files will contain best hits. 
    #uas.txt files will contain unique orthogroup names. Heatmap will be 
    #partitioned by species and by gene and will count raw palindrome counts
def build_heatmap(df):
    species_dict = {"Scer": 1, "Spar": 2, "Smik": 3, "Sbay":4, "Cgla":5, 
                    "Scas":6, "Kwal":7,"Klac":8,"Sklu":9,"Agos":10, "Clus":11,
                    "Dhan":12,"Cgui":13,"Ctro":14,"Calb":15,"Cpar":16,"Lelo":17,
                    "Ylip":18, "Anid":19, "Ncra":20,"Sjap":21,"Soct":22,
                    "Spom": 23}
    df_list = []
    row_dict = {}
    for i in range(1,len(list(species_dict.keys()))):
        row_dict["Species"] = list(species_dict.keys())[i]
        for gene in set(df["Gene"]):
            row_dict[gene] = np.array(len(df.query(f"Gene == '{gene}' and Species == '{list(species_dict.keys())[i]}'")))
        df_list.append(row_dict.copy())
        row_dict.clear()
    
    df = pd.DataFrame(df_list)
    df = df.set_index("Species")
    del df.index.name
    df = df.astype(float)
    
    heat = sns.heatmap(df, yticklabels = list(species_dict.keys())[1:])
    cluster = sns.clustermap(df,yticklabels = list(species_dict.keys())[1:] )
    
    
    heat.figure.savefig("results/ortho_heatmap.png")
    cluster.savefig("results/ortho_cluster.png")
    
    return df
        

if __name__ == "__main__":
    mydata = main()
    df = build_heatmap(mydata)
    
    