import os
import glob
import subprocess
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("--pdb", default="1CRN.pdb", help="Input PDB file")
parser.add_argument("--name", default="Crambin", help="Enzyme Name for output")
args = parser.parse_args()

# --- CONFIGURATION ---
# UPDATE THESE PATHS BEFORE RUNNING ON YOUR LOCAL MACHINE
FOLDX_BINARY = "/path/to/your/foldx.exe" 
windows_base_path = "./output_folder"

PDB_FILE = args.pdb
ENZYME_NAME = args.name

final_output_folder = os.path.join(windows_base_path, ENZYME_NAME)

if not os.path.exists(final_output_folder):
    try:
        os.makedirs(final_output_folder)
        print(f"Created output folder: {final_output_folder}")
    except OSError:
        pass

if not os.path.exists(PDB_FILE):
    print(f"Error: PDB file '{PDB_FILE}' not found.")
    exit()

shutil.copy(PDB_FILE, os.path.join(final_output_folder, PDB_FILE))
print(f"Copied Wild Type to: {final_output_folder}")

list_of_files = glob.glob('**/*.fa', recursive=True)
if not list_of_files:
    print("Error: No .fa files found!")
    exit()

latest_file = max(list_of_files, key=os.path.getctime)
print(f"Processing file: {latest_file}")
print(f"Target Enzyme: {ENZYME_NAME} ({PDB_FILE})")

sequences = []
with open(latest_file, 'r') as f:
    lines = f.readlines()
    for i in range(1, len(lines), 2):
        sequences.append(lines[i].strip())

wild_type_seq = sequences[0]
mutants = sequences[1:]

print(f"Screening {len(mutants)} mutants...")

pdb_base = PDB_FILE.replace(".pdb", "")
foldx_output = f"{pdb_base}_1.pdb"

for i, mutant_seq in enumerate(mutants):
    mutant_id = i + 1
    if len(mutant_seq) != len(wild_type_seq): continue

    mutations_found = []
    is_safe = True
    
    for pos in range(len(wild_type_seq)):
        wt = wild_type_seq[pos]
        mut = mutant_seq[pos]
        if wt != mut:
            code = f"{wt}A{pos + 1}{mut}"
            mutations_found.append(code)
            if mut in ['C', 'M']:
                is_safe = False
                break
    
    if is_safe and mutations_found:
        if os.path.exists(foldx_output):
            os.remove(foldx_output)

        foldx_string = ",".join(mutations_found) + ";"
        print(f"Mutant {mutant_id}: [RUNNING FOLDX]")
        
        with open("individual_list.txt", "w") as f:
            f.write(foldx_string + "\n")
            
        cmd = [
            FOLDX_BINARY,
            "--command=BuildModel",
            f"--pdb={PDB_FILE}",
            "--mutant-file=individual_list.txt"
        ]
        
        subprocess.run(cmd)
        
        if os.path.exists(foldx_output):
            final_name = f"{pdb_base}_{ENZYME_NAME}_Mutant{mutant_id:02d}.pdb"
            destination_path = os.path.join(final_output_folder, final_name)
            
            shutil.move(foldx_output, destination_path)
            print(f"  -> Success! Saved to {destination_path}")
        else:
            print(f"  -> Error: FoldX failed. Check the error message above.")

print("-" * 40)
print("Pipeline Finished.")
