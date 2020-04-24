#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 10:56:01 2020

@author: benjamin
"""


def intersect(files):
    if isinstance(files[0], str):
        
        total_set = {}
        for filename in files:
            with open(filename, "r") as fi:
                for gene in fi.read().split("\n"):
                    if gene in total_set.keys():
                        total_set[gene].append(filename)
                    else:
                        total_set[gene] = []
    return [g for g in total_set.keys() if len(total_set[g]) == len(files)-1]

def union(files):
    total_set = {}
    for filename in files:
        with open(filename, "r") as fi:
            for gene in fi.read().split("\n"):
                total_set[gene] = filename
    return [g for g in total_set.keys()]
                



    
        
            
    