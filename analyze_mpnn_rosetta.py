import pandas as pd
import sys
import os

def analyze_profit(csv_file):
    if not os.path.exists(csv_file):
        print(f"Error: Could not find {csv_file}")
        return

    df = pd.read_csv(csv_file)

    if 'ddG' in df.columns:
        score_col = 'ddG'
    elif 'ddg' in df.columns:
        score_col = 'ddg'
    else:
        print("Error: 'ddG' or 'ddg' column not found. Columns are:")
        print(df.columns.tolist())
        return

    df_sorted = df.sort_values(by=score_col)

    profit_mutants = df_sorted[df_sorted[score_col] < -0.5]
    break_even = df_sorted[(df_sorted[score_col] >= -0.5) & (df_sorted[score_col] < 0)]

    print(f"--- AUDIT REPORT: {csv_file} ---")
    print(f"Total Candidates Screened: {len(df)}")
    print(f"PROFITABLE Mutations ({score_col} < -0.5): {len(profit_mutants)}")
    print(f"Break-even Mutations ({score_col} < 0): {len(break_even)}")

    print("\n--- TOP 5 PROFIT GENERATORS ---")
    if not profit_mutants.empty:
        print(profit_mutants.head(5).to_string(index=False))
        
        max_profit = profit_mutants.iloc[0][score_col]
        budget = abs(max_profit)
        print(f"\nMAXIMUM STABILITY BUDGET: {budget:.2f} REU")
        print("You can afford an Armor mutation that costs up to this amount.")
    else:
        print("None found.")
        print("\nCRITICAL: No profit found. Bankrupt.")

if __name__ == "__main__":
    target_csv = "validation_results.csv"
    if len(sys.argv) > 1:
        target_csv = sys.argv[1]
    
    analyze_profit(target_csv)
