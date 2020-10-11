###############################
# Create module for reading in files from command line

#!/usr/bin/python3
import sys

import Bio.GenBank
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import os 
import sys, getopt
import re

def main():
    input_1=sys.argv[1]
    input_2=sys.argv[2]
    #output=sys.argv[3]
    #print(input_1,input_2,output)
    a =  SeqIO.read(input_1, 'genbank')
    b = SeqIO.read(input_2, "genbank")
    print(a)
    
    
    
    def detect_recombinase_sites(input):
        name=input.seq
        total_out=[]
    

    
        # Recombination sites for the supported recombinases
        phic31={"attP" : "ACGCCCCCAACTGAGAGAACTCG",
                "attB" : "GGAGTACGCGCCCGGGGAGCCC",
                "enzyme_name":"phiC31"}

        bxb1={"attB":"TCGGCCGGCTTGTCGACGACGGCG",
              "attP":"GTCGTGGTTTGTCTGGTCAACCACCGCG",
              "enzyme_name":"Bxb1"}  
    
    
    
    
    
        # List of the recombinases that are currently supported
        recombinase_list=[phic31,bxb1]
        list_out=[]
        for recombinase in recombinase_list:
        
            for key in recombinase:
                for m in re.finditer(str(recombinase[key]),str(name)):
                    central_dinucleotide=str(name[m.end():m.end()+2])
                    info={'recombinase':recombinase['enzyme_name'],'orientation':'forward','position':(m.end()+2),
                          'dinucleotide':central_dinucleotide,'type':key}
                    list_out.append(info)
                    # Searches for a half recombination sites in the forward orientation
            for key in recombinase:
                for m in re.finditer(str((Seq(recombinase[key])).reverse_complement()),str(name)):
                    central_dinucleotide=str(name[m.start()-2:m.start()].reverse_complement())
                    info={'recombinase':recombinase['enzyme_name'],'orientation':'reverse',
                        'position':(m.start()-2),'dinucleotide':central_dinucleotide,'type':key}
                    list_out.append(info)
                    # Serarches for a half recombination site in the reverse orientation 
        total_out.append(list_out)

        return(total_out)

        

    total_sites=detect_recombinase_sites(a)+detect_recombinase_sites(b)

    pair=[]
    outputs=[]
    seen=set()
    for items in total_sites:
        for key in items:
            if key['dinucleotide'] in seen:
                pair.append(key['dinucleotide'])
            else :
                seen.add(key['dinucleotide'])
                seen.add(key['type'])

    for items in total_sites:
        for key in items:
            if key['dinucleotide'] in pair:
                outputs.append(key)
    print(outputs)
    
    
    

if __name__ == "__main__":
    main()



