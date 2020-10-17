#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 08:49:33 2020

@author: jamie
"""
import Recombinase_site_detector as recombine
import sys 
import os
import argparse
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input values', nargs='+', required=True)
    
    args = parser.parse_args()
    if len(args.input) !=2:
        print('Usage: Recombinase_site_detector.py -i <path/to/batch/dir> <common_input>')
        sys.exit()
  
    path=args.input[0]
    os.chdir(path)
    b=args.input[1]
    file_list = [f for f in os.listdir('.') if os.path.isfile(os.path.join('.', f)) and f.endswith('.gb') and f!=b]
    print(file_list)

    for file in file_list:
        a=file
    #print(a)
        output_name='Recombined ' + str(a)
    #print(output_name)
        sys.argv = ['recombinase_detector.py','-i',a,b,'-o',output_name]
        recombine.main(sys.argv)

if __name__ == "__main__":
    main(sys.argv[1:])
