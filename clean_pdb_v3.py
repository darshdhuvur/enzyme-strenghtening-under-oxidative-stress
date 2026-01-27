import sys
import os

KEEP_COFACTORS = {'HEM', 'MN', 'FE', 'ZN', 'MG', 'NAD', 'NDP', 'K'}

def is_valid_line(line):
    if line.startswith("ATOM"):
        return True
    if line.startswith("HETATM"):
        res_name = line[17:20].strip()
        if res_name in KEEP_COFACTORS:
            return True
    return False

def get_sort_key(line):
    chain_id = line[21]
    try:
        res_num = int(line[22:26].strip())
    except ValueError:
        res_num = 9999
    return (chain_id, res_num)

def clean_pdb(input_file):
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}_cleaned.pdb"
    
    print(f"Processing {input_file}...")
    
    with open(input_file, 'r') as f:
        lines = f.readlines()

    valid_lines = [line for line in lines if is_valid_line(line)]
    
    valid_lines.sort(key=get_sort_key)
    
    with open(output_file, 'w') as f:
        current_chain = None
        atom_count = 0
        
        for line in valid_lines:
            chain_id = line[21]
            
            if current_chain and chain_id != current_chain:
                f.write("TER\n")
            
            f.write(line)
            current_chain = chain_id
            atom_count += 1
        
        f.write("TER\n")
        f.write("END\n")

    print(f"Success! Saved to: {output_file}")
    print(f"Total Atoms: {atom_count}")
    print(f"Whitelisted Cofactors Preserved: {KEEP_COFACTORS}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python clean_pdb_v3.py <pdb_file>")
    else:
        clean_pdb(sys.argv[1])
