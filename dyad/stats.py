#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:21:08 2020

@author: benjamin
"""
from math import exp
from scipy.special import ndtr
import scipy.stats as st
from math import factorial as fa
from math import sqrt
from __main__ import averages as av
from statistics import mean as m
from __main__ import unique_palindrome_count as observed_count
from __main__ import total_scanned_sequence_length as total_length
import numpy

#probability of a complementing pair
Ptf = 2*(m(av["C"])*m(av["G"]) + (m(av["A"]))*(m(av["T"])))
#probability of a two-fold palindrome of a half length n with s matching pairs.
def Ps(n,Ptf,s):
    return ((fa(n)*Ptf**s)*(1-Ptf)**(n-s))/((fa(n-s))*fa(s))

#getting the sum of probabilities for a palindrome of half-length between 10-15
#with 0-2 misspairs
stats_ar  = []
for n in range(10,16):
    for s in range(n-2, n):
        #print(s)
        stats_ar.append(2*Ps(n,Ptf,s))
    


total_p = sum(stats_ar)    

_lambda = total_length*total_p
sigma = sqrt(_lambda*(1-total_p))

Zplus = numpy.array((observed_count - 0.5 - _lambda)/sigma)
Zminus = numpy.array((observed_count + 0.5 - _lambda)/sigma)

 

    
    
    
