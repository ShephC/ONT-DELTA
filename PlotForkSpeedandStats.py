import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats

def process_bed_file(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'read_id', 'full_start', 'full_end', 'strand', 'edu_length', 'brdu_length', 'total_length'])
    df['fork_speed'] = df['total_length']/25
    return df

right_forks = process_bed_file('/Users/shephc3/Documents/rightForksFiltered.bed')
left_forks = process_bed_file('/Users/shephc3/Documents/leftForksFiltered.bed')

right_forks['direction'] = 'right'
left_forks['direction'] = 'left'

all_forks = pd.concat([right_forks, left_forks])

stats = all_forks.groupby('direction')['fork_speed'].describe()
print("Fork Speed Statistics:")
print(stats)

t_stat, p_value = scipy_stats.ttest_ind(left_forks['fork_speed'], right_forks['fork_speed'])
print(f"\nt-test results: t-statistic = {t_stat}, p-value = {p_value}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.boxplot([right_forks['fork_speed'], left_forks['fork_speed']])
ax1.set_xticklabels(['Left', 'Right'])
ax1.set_title('Fork Speed Comparison: Left vs Right Progressing Forks')
ax1.set_ylabel('Fork Speed (bp/min)')

ax2.violinplot([left_forks['fork_speed'], right_forks['fork_speed']])
ax2.set_xticks([1, 2])
ax2.set_xticklabels(['Left', 'Right'])
ax2.set_title('Fork Speed Distribution: Left vs Right Progressing Forks')
ax2.set_ylabel('Fork Speed (bp/min)')

plt.tight_layout()
plt.show()