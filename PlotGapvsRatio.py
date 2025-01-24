#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 14:48:42 2025

@author: shephc3
"""

import pandas as pd
import matplotlib.pyplot as plt

def process_bed_files(brdu_file, edu_file):
    # Function to parse a single bed file
    def parse_bed_file(file_path):
        data = {}
        with open(file_path, 'r') as file:
            for line in file:
                fields = line.strip().split()
                read_id = fields[3]
                start = int(fields[1])
                end = int(fields[2])
                if read_id not in data:
                    data[read_id] = []
                data[read_id].append((start, end))
        return data

    # Parse both bed files
    brdu_data = parse_bed_file(brdu_file)
    edu_data = parse_bed_file(edu_file)

    # Collect data for each read
    results = []
    for read_id in set(brdu_data.keys()) | set(edu_data.keys()):
        brdu_tracks = brdu_data.get(read_id, [])
        edu_tracks = edu_data.get(read_id, [])

        if brdu_tracks and edu_tracks:
            brdu_length = sum(end - start for start, end in brdu_tracks)
            edu_length = sum(end - start for start, end in edu_tracks)

            # Calculate distance between tracks
            all_tracks = sorted(brdu_tracks + edu_tracks)
            distances = [all_tracks[i+1][0] - all_tracks[i][1] for i in range(len(all_tracks)-1)]
            min_distance = min(distances) if distances else None

            results.append({
                'Read_ID': read_id,
                'BrdU_Length': brdu_length,
                'EdU_Length': edu_length,
                'Distance': min_distance
            })

    # Create DataFrame
    df = pd.DataFrame(results)
    return df


def plot_gap_vs_ratio(result_df):
    # Calculate the ratio of BrdU to EdU lengths
    result_df['BrdU_EdU_Ratio'] = result_df['BrdU_Length'] / result_df['EdU_Length']

    # Filter out rows with missing or infinite values
    plot_df = result_df.dropna(subset=['Distance', 'BrdU_EdU_Ratio'])
    plot_df = plot_df[plot_df['BrdU_EdU_Ratio'] != float('inf')]

    # Create the scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(plot_df['BrdU_EdU_Ratio'], plot_df['Distance'], alpha=0.5)
    plt.xlabel('BrdU/EdU Length Ratio')
    plt.ylabel('Gap Between Tracks (bp)')
    plt.title('Gap Between Tracks vs BrdU/EdU Length Ratio')
    plt.xscale('log', base=2)  # Use log scale for x-axis due to ratio values
    plt.yscale('linear')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.tight_layout()
    plt.show()

# Usage
brdu_file = '/Users/shephc3/Documents/AllBrdU.bed'
edu_file = '/Users/shephc3/Documents/AllEdU.bed'
result_df = process_bed_files(brdu_file, edu_file)
plot_gap_vs_ratio(result_df)

# Display the first few rows of the DataFrame
print(result_df.head())



# Optionally, save the DataFrame to a CSV file
# result_df.to_csv('processed_data.csv', index=False)
