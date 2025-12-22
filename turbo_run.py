import subprocess
import os
import csv
import glob

foldx_exe = "./foldx_20251231.exe"
pdb_file = "1GV3_clean_Repair.pdb"
mutation_file = "individual_list.txt"
output_csv = "turbo_results.csv"

def run_foldx_turbo():
    if not os.path.exists(pdb_file):
        print(f"CRITICAL ERROR: {pdb_file} not found!")
        return
    if not os.path.exists(mutation_file):
        print(f"Error: {mutation_file} not found. Did you rename it?")
        return

    print(f"--- Starting Turbo Run on {pdb_file} ---")
    print(f"Reading mutations from: {mutation_file}")

    command = [
        foldx_exe,
        "--command=BuildModel",
        f"--pdb={pdb_file}",
        f"--mutant-file={mutation_file}",
        "--numberOfRuns=1", 
        "--pH=7.0"
    ]

    print("Running FoldX... (This might take 10-20 minutes)")
    try:
        subprocess.run(command, check=True)
        print("FoldX run complete!")
    except subprocess.CalledProcessError as e:
        print(f"FoldX crashed! Error code: {e}")
        return

    dif_pattern = f"Dif_{pdb_file.replace('.pdb', '')}*.fxout"
    found_files = glob.glob(dif_pattern)
    
    if not found_files:
        print(f"Error: Could not find any output file matching {dif_pattern}")
        return
    
    dif_filename = found_files[0]
    print(f"Parsing results from {dif_filename}...")
    
    results = []
    with open(dif_filename, 'r') as f:
        lines = f.readlines()
        
        for line in lines:
            if "1GV3_clean_Repair" in line and "total energy" not in line:
                parts = line.split('\t')
                try:
                    mutation_code = parts[1]
                    ddg_value = float(parts[2])
                    results.append([mutation_code, ddg_value])
                except (IndexError, ValueError):
                    continue

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Mutation", "ddG (kcal/mol)"])
        writer.writerows(results)

    print(f"--- SUCCESS ---")
    print(f"Processed {len(results)} mutations.")
    print(f"Results saved to: {output_csv}")

if __name__ == "__main__":
    run_foldx_turbo()
