#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:05:28 2020

@author: benjamin
"""

from fimo_parser.GFF import GFF_Parse
from databases.cadb.CalbDB import CaDB
from databases.scdb.ScerDB import ScDB
from databases.orthdb.OrthDB import OrthDB
from interval_tree import Motif_Node
from dyad_search import UAS
from Bio.Seq import Seq
import re
import time
import pandas as pd




#getting fimo_file
#fimo_file = sys.argv[1]
fimo_file = "wg_fimo/400_fimo_ca_0.01p.gff"
print("extracting fimo hits....")
ts = time.time()
fimo = GFF_Parse(fimo_file)
te= time.time()
print("done in " + str(te-ts) + "s")
df_fimo = fimo.dataframe


#selecting the database 
#database = sys.argv[2]
database = "cadb"

if database == "scdb":
    db = ScDB()
    def _getLocusName(string):
        return string

elif database == "cadb":
    db = CaDB()
    def _getLocusName(string):
        locus_re = re.compile("CAALFM_(.*)")
        code = locus_re.findall(string)[0]
        code_s = code[:2] + "_" + code[2:-1] + "_" + code[-1:len(code)]
        return code_s
elif database == "orthdb":
    db = OrthDB()    



#extracting all coors
s_coor = df_fimo["start"]
e_coor = df_fimo["end"]
locus_re = re.compile("(CAALFM_.*)")
names = [full_name for full_name in df_fimo["name"]]

#getting all sequences
ts = time.time()
print("getting sequences.....")
seq40 = []


#Withdrawing sequence coordinates to scan. If the window of 40bp that I draw 
#around each hit is overlapping with another then I merge the two.
total_scanned_sequence_length = 0
for name in names:
    nodes = [] 
    for row in df_fimo[df_fimo["name"] == name].iterrows():
        #dealing with edge cases
        if row[1]["start"] < 20:
            
            nodes.append(Motif_Node(0, row[1]["end"] + 20))
        elif row[1]["end"] >= 380:
            nodes.append(Motif_Node(row[1]["start"], 400))
        else:
            nodes.append(Motif_Node(row[1]["start"] - 20, row[1]["end"] + 20))
      #decompose() is the function that merges all overlapping features      
    total_promoter_seq = Motif_Node.decompose(nodes)
    total_scanned_sequence_length += sum([node.e - node.s for node in total_promoter_seq])
    for seq in total_promoter_seq:
        seq40.append((name,Seq(db.getSeqFromPromoter(name, 
                       seq.s + 600 , seq.e + 600 ))))   
te = time.time()
print("done in " + str(te-ts) + "s")



print("finding palindromes......")
ts = time.time()
#To find the frequency of each residue
averages = {"C":[], "G":[], "T":[], "A":[]}
unique_palindrome_count = 0
#genes holds all genes with a promoter containing at least one palindrome with
# a decently high cg content to avoid picking up TATA boxes
genes = []
for name, seq in seq40:
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



##Building the output csv and txt documents
#Unique set of gene names
n = set([x.name for x,y in genes])
gene_names = []
pal_dict = {}
csv_out = []
for x,p in genes:
    #print(str(p) + " " + str(p.score))
    
    name = db.lookup(Locus_ID = p.name)[0][0]
    
    if name != "": 
    #print(name)
        pal_dict["Gene"] = name
    else:
        pal_dict["Gene"] = p.name
    
    pal_dict["Sequence"] = p.pal
    pal_dict["Length"] = len(p.pal)
    pal_dict["Score"] = p.score
    csv_out.append(pal_dict.copy())
    pal_dict.clear()
    
csv = pd.DataFrame(csv_out).drop_duplicates()

csv.to_csv("palindromes_out.csv")
#print(len(n))

with open("uas_gene_names.txt", "w") as fo:
    for name in n:
        fo.write(name + "\n")
        

db.close()


        











