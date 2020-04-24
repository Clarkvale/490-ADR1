#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 14:46:17 2020

@author: benjamin
"""

import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bs4 import BeautifulSoup
from databases.orthdb.orthDB import OrthDB
from results.set_operations import union
import time
from Bio.Seq import Seq

class OrthogroupError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None
    def __str__(self):
        if self.message:
            return f"OrthogroupError, {self.message}"
        else:
            return "OrthogroupError raised"

def pull_orthologs(gene):
    url_ortho = """https://portals.broadinstitute.org/cgi-bin/regev/orthogroups/browse_orthogroups.cgi"""
    params = {"orf" : gene, "og_key" : "", "species" : "Any"}
    response = requests.get(url_ortho, params)
    
    assert response.ok, f"{gene} not found"
        
    soup = BeautifulSoup(response.text, "html.parser")
    try:
        og_key =  soup.find_all("font", color = "green")[2].text
        
    except IndexError:
        raise OrthogroupError(f"{gene} not valid query")
    
    
    
    
    sequence_url = f"https://portals.broadinstitute.org/regev/orthogroups/promoters/{og_key[-2]}/{og_key[-1]}/{og_key}-promoters.fasta"
    try:
        seq_request = requests.get(sequence_url)
        
    except:
        raise OrthogroupError(f"Something is wrong with the {gene} query")
    
    assert seq_request.ok, f"{gene} not found"
    
    seqs = seq_request.text.split("\n")
    seq_records = [SeqRecord(seqs[i+1],seqs[i]) for i in range(0,len(seqs) - 1, 2)]

    return seq_records


    
    
if __name__ == "__main__":
    files = ["results/ca_common_name_uas.txt","results/sc_common_name_uas.txt"]
    u =  union(files)
    the_db = OrthDB()
    for name in u:
        try:
            records = pull_orthologs(name)
        except OrthogroupError as e:
            print(e)
            pass
        for record in records:
            the_db.add_seq(name, record)
        time.sleep(10)
     the_db.close()
        
        
        
    the_db = OrthDB()    
    the_db.cursor.execute("SELECT * FROM Promoters")
    sr_list = []
    index = 1
    for row in the_db.cursor.fetchall():
        
        header =  f"{str(index)}_{row[1]}|{row[0]}"
        sr = SeqRecord(seq = Seq(row[2]), id = header)
        sr_list.append(sr)
        index += 1
    with open("fastas/all_valid_orthologs.fasta", "w") as fo:
        SeqIO.write(sr_list, fo, "fasta")
    
    
    the_db.close()
    
    
        
        
    