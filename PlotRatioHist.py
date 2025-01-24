import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def process_bed_file(file_path):
    columns = ['chromosome', 'start', 'end', 'read_id', 'read_start', 'read_end', 'strand',
               'fork_length', 'first_analogue_length', 'second_analogue_length',
               'second_in_first_freq', 'first_in_first_freq',
               'first_in_second_freq', 'second_in_second_freq', 'stall_score']
    
    df = pd.read_csv(file_path, sep='\s+', header=None, names=columns)
    df = df[df['stall_score'] != -3.000000]
    df['ratio'] = df.apply(lambda row: row['second_analogue_length'] / row['first_analogue_length'] 
                           if row['first_analogue_length'] != 0 else np.inf, axis=1)
    return df

left_forks = process_bed_file('/Users/shephc3/Documents/leftForks_stressSignatures.bed')
right_forks = process_bed_file('/Users/shephc3/Documents/rightForks_stressSignatures.bed')

all_forks = pd.concat([left_forks, right_forks])

zero_first_analogue = all_forks[all_forks['first_analogue_length'] == 0]
print(f"Number of reads with zero first analogue length: {len(zero_first_analogue)}")

all_forks = all_forks[all_forks['ratio'] != np.inf]

if all_forks.empty:
    print("No valid data to plot after removing infinite values.")
else:
    plt.figure(figsize=(10, 6))
    plt.hist(all_forks['ratio'], bins=50, edgecolor='black', range=(0, all_forks['ratio'].quantile(0.99)))
    plt.title('Distribution of Second/First Analogue Length Ratios')
    plt.xlabel('Second/First Analogue Ratio')
    plt.ylabel('Frequency')
    plt.axvline(all_forks['ratio'].median(), color='r', linestyle='dashed', linewidth=2)
    plt.savefig('analogue_ratio_distribution.png')
    plt.show()

    print(f"Total forks analyzed: {len(all_forks)}")
    print(f"Mean ratio: {all_forks['ratio'].mean():.2f}")
    print(f"Median ratio: {all_forks['ratio'].median():.2f}")