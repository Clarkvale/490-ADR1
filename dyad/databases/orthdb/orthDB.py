#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 20:51:49 2020

@author: benjamin
"""
from Bio import SeqRecord
import sqlite3

class OrthDB:
    
    def __init__(self, path = "databases/orthdb/Orthologs.db"):
        
        self.connection = sqlite3.connect(path)
        self.cursor = self.connection.cursor()
        self.cursor.execute("""CREATE TABLE IF NOT EXISTS Promoters 
                                (Orthogroup TEXT, Species TEXT, Sequence TEXT)""")
    def add_seq(self, orthogroup, fasta):
        species = fasta.id.split("|")[0][1:]
        self.cursor.execute("""INSERT INTO Promoters VALUES(?,?,?)"""
                            ,(orthogroup,species,fasta.seq))
        self.connection.commit()
        
    def close(self):
        self.cursor.close(); self.connection.close()
      
    def lookup(self,return_q = "*", **kwargs):
        out = []
        for key, value in kwargs.items():
            
            self.cursor.execute(f"SELECT {return_q} FROM Promoters WHERE {key} = ?", 
                                (str(value),))
            out.append(self.cursor.fetchone())
        
        return out 
    
    def lookup_index(self, index):
        self.cursor.execute(f"SELECT * FROM Promoters LIMIT 1 OFFSET {index}")
        return self.cursor.fetchone()
    
    def lookup_seq(self , species, orthogroup):
        out = [] 
        
        self.cursor.execute(f"""SELECT Sequence FROM Promoters WHERE Species = ?
                            AND Orthogroup = ?""", (species, orthogroup))
        out.append(self.cursor.fetchone())
        return out 
    
    def getSeqFromPromoter(self, species, orthogroup, start, end, relative = True):
        query = self.lookup_seq(species, orthogroup)[0]
        
        promoter = query[0]
        if relative == True:
            return promoter[start:end]
        else:
            
            p_sta = query[2]
            rel_start = start - p_sta
            rel_end = (end - start) + rel_start
            return promoter[rel_start:rel_end]
        
if __name__ == "__main__":
    the_db = OrthDB()
    s = the_db.lookup_index("9440")
    
    
