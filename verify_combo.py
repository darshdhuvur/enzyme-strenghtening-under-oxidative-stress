import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
import sys

pyrosetta.init("-beta_nov16 -corrections::beta_nov16 -score:weights beta_nov16_cart -relax:cartesian")

def get_three_letter_code(one_letter):
    mapping = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    return mapping.get(one_letter, 'ALA')

def verify_combo(pdb_file, mutations_list):
    pose = pose_from_pdb(pdb_file)
    scorefxn = create_score_function("beta_nov16_cart")
    
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.cartesian(True)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)

    print(f"--- Building SUPER COMBO: {', '.join(mutations_list)} ---")
    
    for mut in mutations_list:
        old_res = mut[0]
        new_res_char = mut[-1]
        seq_pos = int(mut[1:-1])
        
        new_res_3 = get_three_letter_code(new_res_char)
        
        pose_idx = pose.pdb_info().pdb2pose('A', seq_pos)
        
        if pose_idx == 0:
            print(f"Error: Residue {seq_pos} not found!")
            return None
            
        print(f"Mutating {old_res}{seq_pos} -> {new_res_3} (Pose ID: {pose_idx})")
        pyrosetta.rosetta.protocols.simple_moves.MutateResidue(pose_idx, new_res_3).apply(pose)

    print("Relaxing Combo (this may take ~5-10 mins)...")
    relax.apply(pose)
    
    final_score = scorefxn(pose)
    print(f"Final Combo Score: {final_score:.2f} REU")
    
    mut_str = "_".join(mutations_list)
    filename = f"SuperMutant_{mut_str}.pdb"
    pose.dump_pdb(filename)
    print(f"Saved: {filename}")
    
    return final_score

if __name__ == "__main__":
    INPUT_PDB = "1DGB_cleaned.pdb" 
    
    if len(sys.argv) > 1:
        mutations = sys.argv[1:]
    else:
        mutations = ["H163N", "A21P", "M53L"]

    verify_combo(INPUT_PDB, mutations)
