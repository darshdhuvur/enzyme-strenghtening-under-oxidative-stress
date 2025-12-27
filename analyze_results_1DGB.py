import pandas as pd
import matplotlib.pyplot as plt
import os

INPUT_FILE = "mutations.csv"
OUTPUT_FILE = "top_candidates_1DGB.csv"
PLOT_FILE = "distribution.png"

def analyze():
    if not os.path.exists(INPUT_FILE):
        print(f" Error: Could not find {INPUT_FILE}. Make sure it is in this folder.")
        return

    print(f"reading {INPUT_FILE}...")
    
    try:
        df = pd.read_csv(INPUT_FILE)
    except Exception as e:
        print(f" Error reading CSV: {e}")
        return

    df.columns = df.columns.str.strip()
    if 'ddG' not in df.columns:
        found = False
        for col in df.columns:
            if 'ddg' in col.lower():
                df.rename(columns={col: 'ddG'}, inplace=True)
                found = True
                break
        if not found:
            print(f" Error: Could not find a 'ddG' column in your CSV. Found: {df.columns}")
            return

    df_sorted = df.sort_values(by='ddG', ascending=True)
    
    stabilizing = df_sorted[df_sorted['ddG'] < 0]
    
    print(f"\n Analysis Complete:")
    print(f"   Total Mutations Scanned: {len(df)}")
    print(f"   Stabilizing Mutations (< 0 REU): {len(stabilizing)}")
    print("-" * 40)
    print(" TOP 10 CANDIDATES (Lowest ddG is Best)")
    print("-" * 40)
    print(df_sorted.head(10).to_string(index=False))
    print("-" * 40)

    top_20 = df_sorted.head(20)
    top_20.to_csv(OUTPUT_FILE, index=False)
    print(f" Saved top 20 candidates to '{OUTPUT_FILE}'")

    plt.figure(figsize=(10, 6))
    plt.hist(df['ddG'], bins=50, color='skyblue', edgecolor='black', alpha=0.7)
    plt.axvline(0, color='red', linestyle='dashed', linewidth=1, label='Neutral (0 ΔG)')
    plt.xlabel('ddG (Gibbs Free Energy)')
    plt.ylabel('Count')
    plt.title('Distribution of Mutational Effects (Catalase)')
    plt.legend()
    plt.grid(axis='y', alpha=0.5)
    
    plt.savefig(PLOT_FILE)
    print(f" Saved distribution graph to '{PLOT_FILE}'")
    print("   (Open this image to see the 'Shape of Stability')")

if __name__ == "__main__":
    analyze()
