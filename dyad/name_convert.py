#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:26:22 2020

@author: benjamin
"""
from results.set_operations import intersect, union
from databases.cadb.CalbDB import CaDB
from databases.scdb.ScerDB import ScDB

files = ["results/ca_uas_gene_names.txt","results/sc_uas_gene_names.txt"]

for name in files:
    with open(name, "r") as fi:
        genes = fi.read().split("\n")
    
    if name.split("/")[1].startswith("ca"):
        db = CaDB()
        common_names = []
        for gene in genes:
            n = db.lookup(Locus_ID = gene, return_q = "Gene_ID")[0]
            if n is not None:
                if n[0] != "":
                    common_names.append(n[0])
       
        with open("results/ca_common_name_uas.txt","w") as fo:
            fo.write("\n".join(common_names))
        db.close()
            
    elif name.split("/")[1].startswith("sc"):
        db = ScDB()
        common_names = []
        for gene in genes:
            n = db.lookup(Locus_ID = gene, return_q = "Gene_ID")[0]
            if n != "" and n is not None:
                common_names.append(n[0])
        with open("results/sc_common_name_uas.txt","w") as fo:
            fo.write("\n".join(common_names))
        db.close()
        
