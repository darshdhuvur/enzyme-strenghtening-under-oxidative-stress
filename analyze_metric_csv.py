import pandas as pd
import numpy as np
import sys
import os
from scipy import stats

def calculate_statistics(data):
    sorted_data = np.sort(data)
    n = len(data)
    
    mean = np.mean(data)
    median = np.median(data)
    std = np.std(data, ddof=1)
    sem = stats.sem(data)
    
    q25 = np.percentile(data, 25)
    q75 = np.percentile(data, 75)
    iqr = q75 - q25
    
    iqm_data = data[(data >= q25) & (data <= q75)]
    iqm = np.mean(iqm_data)
    
    trim10 = stats.trim_mean(data, 0.10)
    trim20 = stats.trim_mean(data, 0.20)
    
    filtered_10 = data[(data >= -10) & (data <= 10)]
    filtered_15 = data[(data >= -15) & (data <= 15)]
    mean_filt10 = np.mean(filtered_10) if len(filtered_10) > 0 else np.nan
    mean_filt15 = np.mean(filtered_15) if len(filtered_15) > 0 else np.nan
    
    bottom_3 = np.mean(sorted_data[:3]) if n >= 3 else np.nan
    bottom_5 = np.mean(sorted_data[:5]) if n >= 5 else np.nan
    bottom_10 = np.mean(sorted_data[:10]) if n >= 10 else np.nan
    
    min_val = np.min(data)
    max_val = np.max(data)
    range_val = max_val - min_val
    
    return {
        'N': n,
        'Mean': mean,
        'Median': median,
        'Std': std,
        'SEM': sem,
        'IQM': iqm,
        'Trim10': trim10,
        'Trim20': trim20,
        'Filt_±10': mean_filt10,
        'Filt_±15': mean_filt15,
        'Bottom_3': bottom_3,
        'Bottom_5': bottom_5,
        'Bottom_10': bottom_10,
        'Q25': q25,
        'Q75': q75,
        'IQR': iqr,
        'Min': min_val,
        'Max': max_val,
        'Range': range_val,
        'N_Filt10': len(filtered_10),
        'N_Filt15': len(filtered_15)
    }

def analyze_csv(filepath, mutation_col='Mutation', ddg_col='ddG', 
                wt_score_col='WT_Score', mut_score_col='Mut_Score'):
    
    try:
        df = pd.read_csv(filepath)
        if any(' ' in str(col) for col in df.columns):
            df = pd.read_csv(filepath, sep=r'\s+')
    except Exception as e:
        try:
            df = pd.read_csv(filepath, sep=r'\s+')
        except Exception as e2:
            print(f"Error reading CSV: {e2}")
            return None
    
    df.columns = df.columns.str.strip()
    
    if mutation_col not in df.columns or ddg_col not in df.columns:
        print(f"Error: Required columns not found.")
        print(f"Available columns: {df.columns.tolist()}")
        return None
    
    mutations = df[mutation_col].unique()
    
    results = []
    
    for mut in sorted(mutations):
        mut_data = df[df[mutation_col] == mut]
        ddg_values = mut_data[ddg_col].values
        
        stats_dict = calculate_statistics(ddg_values)
        stats_dict['Mutation'] = mut
        
        if wt_score_col in df.columns:
            wt_scores = mut_data[wt_score_col].values
            stats_dict['WT_Std'] = np.std(wt_scores, ddof=1)
            stats_dict['WT_Range'] = np.max(wt_scores) - np.min(wt_scores)
        
        results.append(stats_dict)
    
    results_df = pd.DataFrame(results)
    
    first_cols = ['Mutation', 'N', 'Mean', 'Median', 'IQM', 'Std', 'SEM']
    other_cols = [col for col in results_df.columns if col not in first_cols]
    results_df = results_df[first_cols + other_cols]
    
    return results_df

def main():
    if len(sys.argv) < 2:
        print("Usage: python ddg_analyzer.py <csv_file>")
        print("\nOptional arguments:")
        print("  --mutation-col <n>")
        print("  --ddg-col <n>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    
    base_name = os.path.splitext(os.path.basename(csv_file))[0]
    output_file = os.path.join(os.getcwd(), f"{base_name}_analyzed.csv")
    
    mutation_col = 'Mutation'
    ddg_col = 'ddG'
    
    for i, arg in enumerate(sys.argv):
        if arg == '--mutation-col' and i+1 < len(sys.argv):
            mutation_col = sys.argv[i+1]
        elif arg == '--ddg-col' and i+1 < len(sys.argv):
            ddg_col = sys.argv[i+1]
    
    results = analyze_csv(csv_file, mutation_col=mutation_col, ddg_col=ddg_col)
    
    if results is None:
        sys.exit(1)
    
    print("\n" + "="*100)
    print("ddG STABILITY ANALYSIS")
    print("="*100)
    print(f"\nInput: {csv_file}")
    print(f"Output: {output_file}")
    print(f"Mutations: {len(results)}")
    print("\n" + "="*100)
    
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.float_format', '{:.3f}'.format)
    print(results.to_string(index=False))
    
    results.to_csv(output_file, index=False, float_format='%.3f')
    
    print("\n" + "="*100)
    print("METRIC DEFINITIONS")
    print("="*100)
    print("""
Mean         : Arithmetic mean of all replicates
Median       : 50th percentile value
IQM          : Interquartile mean (average of 25th-75th percentile)
Trim10       : Mean with top/bottom 10% removed
Trim20       : Mean with top/bottom 20% removed
Filt_±10     : Mean of values within ±10 REU
Filt_±15     : Mean of values within ±15 REU
Bottom_3/5/10: Mean of lowest 3/5/10 values
Q25/Q75      : 25th and 75th percentiles
IQR          : Interquartile range (Q75 - Q25)
WT_Std       : Standard deviation of WT_Score across replicates
WT_Range     : Range of WT_Score values
    """)

if __name__ == "__main__":
    main()
