#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 08:49:33 2020

@author: jamie
"""
import Recombinase_site_detector as recombine
import sys 
import os
os.chdir('/Users/jamie/Documents/Coding/Python/recombinase_simulator/recombinase_simulator')
b='d277_attp-landing-pad-pre-integration.gb'

os.chdir('/Users/jamie/Downloads/benchling_export')
file_list = [f for f in os.listdir('.') if os.path.isfile(os.path.join('.', f)) and f.endswith('.gb')]
print(file_list)
print(b)

for file in file_list:
    a=file
    #print(a)
    output_name='Recombined ' + str(a)
    #print(output_name)
    sys.argv = ['recombinase_detector.py','-i',a,b,'-o',output_name]
    recombine.main(sys.argv)

