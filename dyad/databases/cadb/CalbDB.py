#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:06:54 2020

@author: benjamin
"""

import sqlite3

class CaDB:
    
    def __init__(self):
        
        self.connection = sqlite3.connect("databases/cadb/C.albicans_CDS.db")
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

    def getSeqFromPromoter(self, name, start, end):
        query = self.lookup(Locus_ID = name)[0]
    
        promoter = query[4]
        
        p_sta = query[2]
        rel_start = start - p_sta
        rel_end = (end - start) + rel_start
        return promoter[rel_start:rel_end]
    
     
if __name__ == "__main__":
    the_db = CaDB()
    mvd = the_db.lookup(Gene_ID = "MVD")[0]
    seqat200 = the_db.getSeqFromPromoter(mvd[1],12900,12950)
    the_db.close()
    
    
        
        
        