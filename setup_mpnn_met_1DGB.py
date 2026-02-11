import Bio.PDB
import json
import sys
import os

INPUT_TETRAMER_PDB = "1DGB_cleaned.pdb" 
CHAIN_ID = "A"
INTERFACE_CUTOFF = 4.0
MPNN_NAME = "1DGB_cleaned_monomer"
OUTPUT_SEQS = "seqs.jsonl"
OUTPUT_FIXED = "fixed_positions.jsonl"
CATALYTIC_RESIDUES_PDB = [75, 148, 358, 354, 74]

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def get_interface_residues(structure, target_chain_id, cutoff):
    model = structure[0]
    target_chain = model[target_chain_id]
    other_atoms = []
    
    for chain in model:
        if chain.id != target_chain_id:
            for residue in chain:
                if residue.id[0] == ' ':
                    for atom in residue:
                        other_atoms.append(atom)

    if not other_atoms:
        print("!! WARNING: No other chains found. Cannot calculate interface.")
        return set()

    ns = Bio.PDB.NeighborSearch(other_atoms)
    interface_residues = set()
    
    for residue in target_chain:
        if residue.id[0] != ' ': continue
        
        is_interface = False
        for atom in residue:
            if ns.search(atom.coord, cutoff):
                is_interface = True
                break
        
        if is_interface:
            interface_residues.add(residue.id[1])
            
    return interface_residues

def main():
    parser = Bio.PDB.PDBParser(QUIET=True)
    
    if not os.path.exists(INPUT_TETRAMER_PDB):
        print(f"Error: Could not find {INPUT_TETRAMER_PDB}")
        sys.exit(1)

    try:
        structure = parser.get_structure("1DGB", INPUT_TETRAMER_PDB)
    except Exception as e:
        print(f"Error reading structure: {e}")
        sys.exit(1)
    
    print(f"Reading {INPUT_TETRAMER_PDB} to identify interface residues...")
    
    interface_pdb_ids = get_interface_residues(structure, CHAIN_ID, INTERFACE_CUTOFF)
    print(f"Found {len(interface_pdb_ids)} interface residues on Chain {CHAIN_ID}.")

    chain = structure[0][CHAIN_ID]
    
    sequence = ""
    fixed_list = []
    designable_mets = []
    mpnn_index = 1 
    
    print("\nScanning Chain A...")
    print(f"{'PDB ID':<10} | {'Res':<5} | {'MPNN Idx':<10} | {'Status'}")
    print("-" * 50)

    for residue in chain:
        if residue.id[0] != ' ':
            continue
            
        resname_3 = residue.get_resname()
        resname = THREE_TO_ONE.get(resname_3, 'X')
            
        sequence += resname
        pdb_num = residue.id[1]
        
        is_met = (resname == 'M')
        is_interface = (pdb_num in interface_pdb_ids)
        is_catalytic = (pdb_num in CATALYTIC_RESIDUES_PDB)
        
        should_design = is_met and not is_interface and not is_catalytic
        
        status = ""
        if should_design:
            designable_mets.append(f"M{pdb_num}")
            status = "DESIGN (Met)"
        else:
            fixed_list.append(mpnn_index)
            if is_catalytic: status = "FIXED (Catalytic)"
            elif is_interface: status = "FIXED (Interface)"
            else: status = "FIXED (Non-Met)"

        if is_met or is_catalytic:
            print(f"{pdb_num:<10} | {resname:<5} | {mpnn_index:<10} | {status}")

        mpnn_index += 1

    seq_record = {"name": MPNN_NAME, "seq": sequence}
    with open(OUTPUT_SEQS, 'w') as f:
        f.write(json.dumps(seq_record) + "\n")
        
    fixed_record = {MPNN_NAME: {CHAIN_ID: fixed_list}}
    with open(OUTPUT_FIXED, 'w') as f:
        f.write(json.dumps(fixed_record) + "\n")

    print("-" * 50)
    print(f"Total Residues: {len(sequence)}")
    print(f"Fixed Residues: {len(fixed_list)}")
    print(f"Designable Methionines: {len(designable_mets)}")
    print(f"Targets: {', '.join(designable_mets)}")
    print("-" * 50)
    print(f"Saved: {OUTPUT_SEQS}")
    print(f"Saved: {OUTPUT_FIXED}")

if __name__ == "__main__":
    main()
