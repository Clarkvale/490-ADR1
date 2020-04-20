#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 16:11:20 2020

@author: benjamin
"""

from functools import total_ordering 
    
@total_ordering 
class Motif_Node:
    def __init__(self, start, end):
        self.s = start
        self.e = end
     
    def __contains__(self, other_int):
        return ((self.s <= other_int.s and other_int.s <= self.e) 
                or (self.s >= other_int.s and other_int.e >= self.s))
    
    def __add__(self, other):
        new_start = min([self.s,other.s])
        new_end = max([self.e,other.e])
        return Motif_Node(new_start, new_end)
    
    def __lt__(self, other):
        return (self.s < other.s)
    
    def __eq__(self, other):
        return (self.s == other.s)
    
    def __nt__(self, other):
        return not (self.s == other.s)
    def __repr__(self):
        return str((self.s,self.e))
    
    @staticmethod
    def decompose(list_of_nodes):
        sorted_nodes = sorted(list_of_nodes)
        i = 0
        while i < len(sorted_nodes) - 1:
            if sorted_nodes[i] in sorted_nodes[i + 1]:
                sorted_nodes[i] = sorted_nodes[i] + sorted_nodes[i + 1]
                sorted_nodes.pop(i + 1)
            else:
                i += 1
        return sorted_nodes

if __name__ == "__main__":
    nodes = [(50,107),(20,67),(100,147),(220,267),(200,247), (10,19),(10,19)]
    motif_nodes = [Motif_Node(node[0],node[1]) for node in nodes]
    print(Motif_Node.decompose(motif_nodes))
    
    
            
                            
                    
                    
                    

            
        