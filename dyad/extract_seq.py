#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:40:29 2020

@author: benjamin
"""

from databases.cadb.CalbDB import CaDB
from databases.scdb.ScerDB import ScDB
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
##add command-line functionality
##probably just make one function

#candida
#cadb = CaDB()
#cadb.cursor.execute("SELECT * FROM Promoters")
#
#def _getLocusName(string):
#    locus_re = re.compile("CAALFM_(.*)")
#    code = locus_re.findall(string)[0]
#    code_s = code[:2] + "_" + code[2:-1] + "_" + code[-1:len(code)]
#    return code_s
#
#out = []
#for row in cadb.cursor.fetchall():
#    name = _getLocusName(row[1])
#    seq = Seq(row[4][599:999])
#    out.append(SeqRecord(seq, name))
#
#cadb.close()
#
#with open("fastas/ca_400bp_promoter.fasta", "w") as fo:
#    SeqIO.write(out, fo, "fasta")
    

#cerevisiae
scdb  = ScDB()
scdb.cursor.execute("SELECT * FROM Promoters")

out = []
for row in scdb.cursor.fetchall():
    
    name = row[1]
        
    seq = Seq(row[5][600:999])
    out.append(SeqRecord(seq, name))
    
scdb.close()

with open("fastas/sc_400bp_promoter.fasta", "w") as fo:
    SeqIO.write(out, fo, "fasta")

    

