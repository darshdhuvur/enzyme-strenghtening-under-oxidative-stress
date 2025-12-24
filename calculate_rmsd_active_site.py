import Bio.PDB
import numpy as np
import sys
import os

WT_PDB = "1GV3_clean_Repair.pdb"
MUTANT_PDB = "1GV3_triple_mutant_KA52R+GA159L+GA27A_Repair.pdb"

def get_active_site_residues(structure):
    model = structure[0]
    
    mn_atom = None
    active_site_residues = []

    print("  Scanning for Manganese...")
    for chain in model:
        for residue in chain:
            if residue.resname.strip() == "MN":
                for atom in residue:
                    mn_atom = atom
                    print(f"  -> Found MN in Chain '{chain.id}' Residue {residue.id[1]}")
                    break
            if mn_atom: break
        if mn_atom: break
    
    if not mn_atom:
        print("  [Error] No Manganese (MN) atom found in structure!")
        return []

    print("  Scanning for Ligands (< 2.5 A)...")
    for chain in model:
        for residue in chain:
            if residue.resname == "MN": continue
            if residue.resname == "HOH": continue

            for atom in residue:
                diff = atom.coord - mn_atom.coord
                dist = np.sqrt(np.sum(diff * diff))
                
                if dist < 2.5:
                    if residue not in active_site_residues:
                        active_site_residues.append(residue)
                        print(f"    -> Ligand: {residue.resname} {residue.id[1]} (Chain {chain.id})")
                    break
    
    return active_site_residues

def calculate_active_site_rmsd(wt_file, mut_file):
    parser = Bio.PDB.PDBParser(QUIET=True)

    if not os.path.exists(wt_file):
        print(f"Error: Wild Type file not found: {wt_file}")
        return
    if not os.path.exists(mut_file):
        print(f"Error: Mutant file not found: {mut_file}")
        return

    print(f"Loading WT: {wt_file}")
    wt_struct = parser.get_structure("WT", wt_file)
    
    print(f"Loading Mutant: {mut_file}")
    mut_struct = parser.get_structure("Mutant", mut_file)

    print("\n--- Identifying Active Site (WT) ---")
    wt_residues = get_active_site_residues(wt_struct)
    
    print("\n--- Identifying Active Site (Mutant) ---")
    mut_residues = get_active_site_residues(mut_struct)

    if len(wt_residues) == 0:
        print("Error: Could not find active site in WT.")
        return

    wt_atoms = []
    mut_atoms = []

    for wt_res in wt_residues:
        res_num = wt_res.id[1]
        chain_id = wt_res.get_parent().id
        
        match = None
        for mut_res in mut_residues:
            if mut_res.id[1] == res_num and mut_res.get_parent().id == chain_id:
                match = mut_res
                break
        
        if match and 'CA' in wt_res and 'CA' in match:
            wt_atoms.append(wt_res['CA'])
            mut_atoms.append(match['CA'])
        else:
            print(f"Warning: Could not match CA atoms for residue {res_num}")

    if len(wt_atoms) == 0:
        print("Error: No common atoms found for RMSD calculation.")
        return

    superimposer = Bio.PDB.Superimposer()
    superimposer.set_atoms(wt_atoms, mut_atoms)
    
    rmsd = superimposer.rms
    print(f"\nRESULT Active Site RMSD: {rmsd:.5f} Angstroms")

    if rmsd <= 0.2:
        print("STATUS: PRESERVED (Ideal)")
    elif rmsd <= 0.6:
        print("STATUS: TOLERABLE (Caution)")
    else:
        print("STATUS: DISTORTED (Fail)")

if __name__ == "__main__":
    calculate_active_site_rmsd(WT_PDB, MUTANT_PDB)
