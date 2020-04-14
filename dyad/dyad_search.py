#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 14:04:40 2020

@author: benjamin
"""


from Bio.Seq import Seq


        
#Upstream Activation Sequence
class UAS:
    
    class _Palindrome:
        def __init__(self, name, pal, score, pos):
            self.name = name
            self.pal = pal
            self.score = score
            self.pos = pos 
            
        def __repr__(self):
            return str(self.pal)
        
        def __str__(self):
            return str(self.pal)
        
        def cg_content(self):
            return ((self.pal.count("C") + self.pal.count("G")) /len(self.pal))*100
        
        #if these are confusing its because they are. Im using the heuristic that 
        #"greater" palindromes are better, which means they are both larger 
        #and have a lower score. The score always takes precedent, meaning
        #the size is used as a tie breaker if the score is the same. 
        
        def __len__(self):
            return len(self.pal)
        
        def __eq__(self, other):
            return len(self.pal) == len(other.pal) and self.score == other.score 
        
        def __ne__(self, other):
            return self.pal != other.pal or self.score != other.score
        
        def __lt__(self, other):
            if self.score == other.score:
                return len(self.pal) < len(other.pal)
            else:
                return self.score > other.score
            
        def __le__(self, other):
            if self.score == other.score:
                return len(self.pal) <= len(other.pal)
            else:
                return self.score >= other.score
            
        def __gt__(self, other):
            if self.score == other.score:
                return len(self.pal) > len(other.pal)
            else:
                return self.score < other.score
            
        def __ge__(self,other):
            if self.score == other.score:
                return len(self.pal) >= len(other.pal)
            else:
                return self.score <= other.score
                
                
                
    
    def __init__(self, name, string, missmatches = 3, size_range = (20,30)):
        self.name = name
        self.min_size, self.max_size = size_range
        if not isinstance(string,Seq):
            self.sequence = Seq(string.upper())
        else:
            self.sequence = string
        self.pals_in_line = []
        
        self.__find_pals(missmatches)
        
    def __repr__(self):
        return str((self.name, self.pals_in_line))
        

        
    def __find_pals(self, n):
        for k in range(self.min_size, self.max_size, 2):
            for kmer, pos in self._window_slide(self.sequence,k):
                ispal  = self.is_palindrome(kmer, n)
                if ispal[0]:
                    self.pals_in_line.append(self._Palindrome(self.name,kmer,ispal[1], pos))
                    
    def get_best_pal(self):
       return max(sorted(self.pals_in_line))
        
    def is_palindrome(self, string, n):
        assert len(string) % 2 == 0, "Sequence not divisible by 2."
        COMPLEMENT = {"G":"C", "C":"G", "T":"A", "A":"T"}
        m = 0
        if not isinstance(string, Seq):
            seq = Seq(string)
        else:
            seq = string
            
        if seq == seq.reverse_complement():
            return (True,m)
        
           
        for i in range(int(len(seq)/2)):
            try:
                if seq[i] != COMPLEMENT[seq[len(seq)- i -1]]:
                    m += 1
                    if m > n:
                        return (False,m)
            except KeyError:
                m += 1
                if m > n:
                        return (False,m)
        return (True,m)

    def _window_slide(self, seq, k):
        for i in range(len(seq)-k + 1):
            t = (seq[i: i + k], i)
            yield t


        
if __name__ == "__main__":
    
    seq = "TTTGTCTTCATTAACGGCTTTCGCTCATAAAAATGTTATGACGTTTTGCCCGCAGGCGGGAAACCATC"
    + "CACTTCACGAGACTGATCTCCTCTGCCGGAACACCGGGCATCTCCAACTTATAAGTTGGAGAAATAAGAGAA"
    + "TTTCAGATTGAGAGAATGAAAAAAAAAAAAAAAAAAAAGGCAGAGGAGAGCATAGAAATGGGGTTCACTTTT"
    + "TGGTAAAGCTATAGCATGCCTATCACATATAAATAGAGTGCCAGTAGCGACTTTTTTCACACTCGAAATACT"
    + "CTTACTACTGCTCTCTTGTTGTTTTTATCACTTCTTGTTTCTTCTTGGTAAATAGAATATCAAGCTACAAAA"
    + "AGCATACAATCAACTATCAACTATTAACTATATCGTAATACACA"
    
    uas = UAS("dyad", seq, 0)
    print(uas.get_best_pal())
    print(len(uas.get_best_pal()))
#    locus_re = re.compile("(CAALFM_.*)")
#    
#    #dyad search
#    erg_fimo = GFF_Parse("wg_fimo/fimo.gff")
#    erg_df = erg_fimo.dataframe
#    
#    end_coor = erg_df["end"]
#    start_coor = erg_df["start"]
#    
#    names = [locus_re.findall(full_name)[0] for full_name in erg_df["name"]]
#    
#    seq40= []
#    db = CaDB()
#    for i in range(len(names)):
#        seq40.append((names[i],Seq(db.getSeqFromPromoter(names[i], 
#                       start_coor[i] - 15  , end_coor[i] + 15 ))))
#    db.close()
#    
#    uas1 = UAS(seq40[0][0], seq40[0][1])
#    
#    genes = []
#    for name, seq in seq40:
#        uas = UAS(name, seq, 3)
#        if any(uas.pals_in_line):
#            genes.append((uas, uas.get_best_pal()))
#    
#    
#    #cg content filter
#    cg_genes = []
#    for uas, pal in genes:
#        if uas.cg_content() > 20 and uas.cg_content() < 35:
#            cg_genes.append((uas,pal))
            
    #
            
        
    
            
    
        

    
     
     
     
 




