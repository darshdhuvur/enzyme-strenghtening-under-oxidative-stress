import sys

def parse_pdb_residues(pdb_file):
    weak_spots = []
    aa_map = {'CYS': 'C', 'MET': 'M'}

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and "CA" in line:
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = line[22:26].strip()
                
                if res_name in aa_map:
                    code = aa_map[res_name]
                    weak_spots.append((code, chain, res_num))
    return weak_spots

def generate_foldx_list(weak_spots, output_file="individual_list.txt"):
    with open(output_file, 'w') as f:
        all_mutations = []
        for code, chain, res_num in weak_spots:
            new_aa = 'S' if code == 'C' else 'L'
            
            mutation_str = f"{code}{chain}{res_num}{new_aa}"
            all_mutations.append(mutation_str)
            
            f.write(f"{mutation_str};\n")
        
        if all_mutations:
            f.write(f"{','.join(all_mutations)};\n")

    print(f"Found {len(weak_spots)} weak spots.")
    print(f"Saved mutation list to: {output_file}")

if __name__ == "__main__":
    pdb_filename = "1ubq.pdb"
    
    try:
        print(f"Scanning {pdb_filename}...")
        targets = parse_pdb_residues(pdb_filename)
        generate_foldx_list(targets)
    except FileNotFoundError:
        print(f"Error: Could not find {pdb_filename}")
