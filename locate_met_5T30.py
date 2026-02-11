import Bio.PDB
import sys
import os

INPUT_PDB = "5T30_cleaned.pdb"
CHAIN_ID = "A"
INTERFACE_CUTOFF = 6.0

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
        return set()

    ns = Bio.PDB.NeighborSearch(other_atoms)
    interface_residues = set()
    
    for residue in target_chain:
        if residue.id[0] != ' ': continue
        
        for atom in residue:
            if ns.search(atom.coord, cutoff):
                interface_residues.add(residue.id[1])
                break
            
    return interface_residues

def main():
    parser = Bio.PDB.PDBParser(QUIET=True)
    
    if not os.path.exists(INPUT_PDB):
        print(f"Error: Could not find {INPUT_PDB}")
        sys.exit(1)

    try:
        structure = parser.get_structure("5T30", INPUT_PDB)
    except Exception as e:
        print(f"Error reading structure: {e}")
        sys.exit(1)
    
    interface_pdb_ids = get_interface_residues(structure, CHAIN_ID, INTERFACE_CUTOFF)

    chain = structure[0][CHAIN_ID]
    designable_mets = []

    for residue in chain:
        if residue.id[0] != ' ':
            continue
            
        resname = residue.get_resname()
        pdb_num = residue.id[1]
        
        if resname == 'MET':
            if pdb_num not in interface_pdb_ids:
                designable_mets.append(f"M{pdb_num}")

    print(f"Found {len(designable_mets)} mutatable methionines (Non-Interface > {INTERFACE_CUTOFF}A):")
    print(", ".join(designable_mets))

if __name__ == "__main__":
    main()
