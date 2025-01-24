import pandas as pd

def find_overlaps(group):
    overlaps = []
    for i in range(len(group) - 1):
        current = group.iloc[i]
        next_read = group.iloc[i + 1]
        if current['Track End'] >= next_read['Track Start']:
            overlaps.append((
                str(current['Chromosome']),
                int(current['Track Start']),
                int(current['Track End']),
                str(current['Read ID']),
                int(current['Read Start']),
                int(current['Read End']),
                str(current['Strand'])
            ))
            overlaps.append((
                str(next_read['Chromosome']),
                int(next_read['Track Start']),
                int(next_read['Track End']),
                str(next_read['Read ID']),
                int(next_read['Read Start']),
                int(next_read['Read End']),
                str(next_read['Strand'])
            ))
    return overlaps

columnNames = ['Chromosome', 'Track Start', 'Track End', 'Read ID', 'Read Start', 'Read End', 'Strand']
originDF = pd.read_csv("/Users/shephc3/Documents/origins_unique.bed", sep='\t', names=columnNames)

df_sorted = originDF.sort_values(['Chromosome', 'Track Start'])

overlapping_reads = df_sorted.groupby('Chromosome').apply(find_overlaps)
all_overlaps = [item for sublist in overlapping_reads for item in sublist]

output_df = pd.DataFrame(all_overlaps, columns=columnNames)
output_df = output_df.drop_duplicates()

output_df.to_csv("/Users/shephc3/Documents/overlapping_reads.bed", sep='\t', header=False, index=False)

print(f"BED file exported with {len(output_df)} overlapping reads.")