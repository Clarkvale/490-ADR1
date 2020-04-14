#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 14:38:50 2020

@author: benjamin
"""
import pandas as pd
import re

class GFF_Parse:
    
    def __tokenize(self, line):
        line_list = line.split("\t")
        

        line_dict = {}
        
        line_dict["header"] = line_list[0]
        
        line_dict["start"] = int(line_list[3])
        
        line_dict["end"] = int(line_list[4])
        
        line_dict["strand"] = line_list[6]
        
        line_dict["sequence"] = self.rx_dict["sequence"].findall(line_list[8])[0]
        
        line_dict["motif_id"] = self.rx_dict["motif_id"].findall(line_list[8])[0]
        
        line_dict["name"] = self.rx_dict["name"].findall(line_list[8])[0]
        
        try:
           
            line_dict["pvalue"] = float(self.rx_dict["pvalue"].findall(line_list[8])[0])
            
        except IndexError:
            
            if "pvalue" not in line_dict.keys():
                line_dict["pvalue"] = float(self.rx_dict["psci_val"].findall(line_list[8])[0])
            
        try:
             line_dict["qvalue"] = float(self.rx_dict["qvalue"].findall(line_list[8])[0])
             
        except IndexError:
             
             if "qvalue" not in line_dict.keys():
                line_dict["qvalue"] = float(self.rx_dict["qsci_val"].findall(line_list[8])[0])
                
            
        
        return line_dict
    
    
    def __parse(self, pathToFile):
        data_list = []
        with open(pathToFile, 'r') as fi:
            line = fi.readline()
            
            while line:
                
                if not line.startswith("##"):
                    tokens = self.__tokenize(line)
                    data_list.append(tokens)
                    line = fi.readline()
                else:
                    line = fi.readline()
        
                
            
        return data_list
        
    def __init__(self, pathToFile):
        self.path = pathToFile
        self.rx_dict = {"pvalue" : re.compile(";pvalue=(0.\d+);"),
                    "qvalue": re.compile(";qvalue=(0.\d+);"),
                    "sequence": re.compile(";sequence=(\w+);"),
                    "psci_val" : re.compile(";pvalue=(\d\.?\d*e-\d+);"),
                    "motif_id": re.compile("ID=([^-]+)"),
                    "qsci_val" : re.compile(";qvalue=(\d\.?\d*e-\d+);"),
                    "name" : re.compile("Name=[^_]*_([^;+]*(?<=[^-+]))")}
        

        self.dataframe = pd.DataFrame(self.__parse(self.path))
        self.dataframe = self.dataframe.astype({"header": str})
        
        
if __name__ == "__main__":
    fimo_obj = GFF_Parse("/home/benjamin/repo/common/Bioinfo/490/ADR1/Candida/Erg_pathway/fimo_out_0.001/fimo.gff")
    mydf = fimo_obj.dataframe
    #print(mydf)
       
    
    
    
