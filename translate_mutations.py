import os

pdb_file = "1GV3_MASTER.pdb"
fasta_file = "output_1GV3/seqs/1GV3_MASTER.fa"

aa_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

pdb_residues = []
last_res_num = None
last_chain = None

print("Mapping PDB structure...")
with open(pdb_file, "r") as f:
    for line in f:
        if line.startswith("ATOM") and "CA" in line:
            chain = line[21]
            res_num = line[22:26].strip()
            res_name = line[17:20].strip()
            
            if res_num != last_res_num or chain != last_chain:
                code = aa_map.get(res_name, 'X')
                pdb_residues.append((chain, res_num, code))
                last_res_num = res_num
                last_chain = chain

print(f"Mapped {len(pdb_residues)} residues from PDB.")

sequences = []
with open(fasta_file, "r") as f:
    current_seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq: sequences.append(current_seq)
            current_seq = ""
        else:
            current_seq += line
    if current_seq: sequences.append(current_seq)

wt_seq = sequences[0]
ai_seq = sequences[1] 

final_list = []

print("\nTranslating mutations...")
for i, (wt_aa, ai_aa) in enumerate(zip(wt_seq, ai_seq)):
    if wt_aa != ai_aa:
        if i < len(pdb_residues):
            chain, res_num, pdb_aa = pdb_residues[i]
            
            if wt_aa != pdb_aa:
                print(f"Warning: Mismatch at index {i}. AI thinks it's {wt_aa}, PDB says {pdb_aa}")
                continue
                
            mutation_code = f"{wt_aa}{chain}{res_num}{ai_aa};"
            final_list.append(mutation_code)

print(f"\nGenerated {len(final_list)} valid FoldX mutations.")
print("Saving to mpnn_mutations.txt...")

with open("mpnn_mutations.txt", "w") as f:
    for mut in final_list:
        f.write(mut + "\n")

print("Done! Check mpnn_mutations.txt")
