import os
import glob

list_of_files = glob.glob('**/*.fa', recursive=True)

if not list_of_files:
    print("Error: No .fa files found! Checked current folder and subfolders.")
    exit()

latest_file = max(list_of_files, key=os.path.getctime)
fasta_file = latest_file

print(f"Auto-detected latest ProteinMPNN output: {fasta_file}")
print("-" * 40)

sequences = []
with open(fasta_file, 'r') as f:
    lines = f.readlines()
    for i in range(1, len(lines), 2):
        seq = lines[i].strip()
        sequences.append(seq)

wild_type_seq = sequences[0]
mutants = sequences[1:]

print(f"Found Wild Type (Length: {len(wild_type_seq)})")
print(f"Found {len(mutants)} mutants to screen.")
print("-" * 40)

valid_mutants_for_foldx = []

for i, mutant_seq in enumerate(mutants):
    mutant_id = i + 1
    
    if len(mutant_seq) != len(wild_type_seq):
        continue

    mutations_found = []
    is_radiation_safe = True
    
    for position in range(len(wild_type_seq)):
        wt_aa = wild_type_seq[position]
        mut_aa = mutant_seq[position]
        
        if wt_aa != mut_aa:
            mutation_code = f"{wt_aa}A{position + 1}{mut_aa}"
            mutations_found.append(mutation_code)
            
            if mut_aa == 'C' or mut_aa == 'M':
                is_radiation_safe = False
                break 
    
    if is_radiation_safe and len(mutations_found) > 0:
        foldx_string = ",".join(mutations_found)
        print(f"Mutant {mutant_id}: [ACCEPTED] - {foldx_string}")
        valid_mutants_for_foldx.append(foldx_string)
            
print("-" * 40)
print(f"Screening Complete. {len(valid_mutants_for_foldx)} mutants passed the filter.")
