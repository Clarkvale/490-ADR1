#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 13:44:18 2020

@author: benjamin
"""

import sqlite3

class ScDB:
    
    def __init__(self):
        
        self.connection = sqlite3.connect("databases/scdb/S.cerevisiae_CDS.db")
        self.cursor = self.connection.cursor()
        
    def close(self):
        self.cursor.close(); self.connection.close()
      
    def lookup(self, **kwargs):
        out = []
        for key, value in kwargs.items():
            
            self.cursor.execute(f"SELECT * FROM Promoters WHERE {key} = ?", 
                                (str(value),))
            out.append(self.cursor.fetchone())
        
        return out 
    
    def lookup_seq(self , **kwargs):
        out = [] 
        for key, value in kwargs.items():
            self.cursor.execute(f"SELECT Sequence FROM Promoters WHERE {key} = ?",
                                (str(value),))
            out.append(self.cursor.fetchone())
        return out 
    
    def getSeqFromPromoter(self, name, start, end, relative = True):
        query = self.lookup_seq(Locus_ID = name)[0]
        
        promoter = query[0]
        if relative == True:
            return promoter[start:end]
        else:
            
            p_sta = query[2]
            rel_start = start - p_sta
            rel_end = (end - start) + rel_start
            return promoter[rel_start:rel_end]
    
     
if __name__ == "__main__":
    the_db = ScDB()
    upc2 = the_db.lookup(Gene_ID = "ADH2")[0]
    #print(upc2)
    seqat200 = the_db.getSeqFromPromoter(upc2[1],600,)
    the_db.close()
    
    
        
        
        