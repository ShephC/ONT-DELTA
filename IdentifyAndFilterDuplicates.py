#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 09:38:17 2025

@author: shephc3
"""

import pandas as pd

def check_duplicates(file_path):
    # Read only the first 4 columns of the BED file
    df = pd.read_csv(file_path, sep='\s+', header=None, 
                     usecols=[0, 1, 2, 3, 4, 5, 6],
                     names=['chr', 'start', 'end', 'read_id', 'read_start', 'read_stop', 'strand'])
    
    # Check for duplicate read_ids
    duplicates = df[df.duplicated(subset='read_id', keep=False)]
    
    if duplicates.empty:
        print(f"No duplicate reads found in {file_path}")
    else:
        print(f"Duplicate reads found in {file_path}:")
        print(duplicates)
        print(f"Total number of duplicates: {len(duplicates)}")

# Paths to your BED files
bed_files = [
    '/Users/shephc3/Documents/AllBrdU.bed',
    '/Users/shephc3/Documents/AllEdU.bed',
    '/Users/shephc3/Documents/leftForks_stressSignatures.bed',
    '/Users/shephc3/Documents/rightForks_stressSignatures.bed',
    '/Users/shephc3/Documents/origins.bed',
    '/Users/shephc3/Documents/terminations.bed'
]

# Check each file for duplicates
for file in bed_files:
    check_duplicates(file)


def remove_duplicates_and_write(input_file, output_file):
    df = pd.read_csv(input_file, sep='\s+', header=None, 
                     usecols=[0, 1, 2, 3, 4, 5, 6],
                     names=['chr', 'start', 'end', 'read_id', 'read_start', 'read_stop', 'strand'])
    
    df_unique = df.drop_duplicates(subset='read_id', keep='first')
    
    df_unique.to_csv(output_file, sep='\t', header=False, index=False)
    
    print(f"Processed {input_file}")
    print(f"Original entries: {len(df)}")
    print(f"Unique entries: {len(df_unique)}")
    print(f"Removed {len(df) - len(df_unique)} duplicates")
    print(f"Written to {output_file}\n")


for input_file in bed_files:
    output_file = input_file.replace('.bed', '_unique.bed')
    remove_duplicates_and_write(input_file, output_file)