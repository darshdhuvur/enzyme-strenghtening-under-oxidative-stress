import pyrosetta
from pyrosetta import init, pose_from_pdb, create_score_function
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.core.scoring import rms_at_corresponding_atoms
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID
from pyrosetta.rosetta.utility import vector1_unsigned_long
import multiprocessing
import os
import pandas as pd
import numpy as np
import Bio.PDB
import shutil

PDB_FILE = "1DGB_cleaned_monomer.pdb"
OUTPUT_CSV = "final_1DGB_mutants_results.csv"
THREADS_TO_USE = 12
N_REPLICATES = 25
ACTIVE_SITE_CUTOFF = 5.0

COMBOS = {
    "Quadruple_All": [
        ("M", 287, "L"), 
        ("M", 353, "L"), 
        ("M", 339, "L"), 
        ("T", 43, "R")
    ],
    "Budget_Control": [
        ("M", 339, "L"), 
        ("T", 43, "R")
    ]
}

def init_rosetta():
    try:
        init(
            "-relax:cartesian "
            "-score:weights ref2015_cart "
            "-relax:min_type lbfgs_armijo_nonmonotone "
            "-relax:default_repeats 5 "
            "-ex1 -ex2 "
            "-mute all"
        )
    except RuntimeError:
        pass

def get_3_letter_code(one_letter):
    mapping = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    return mapping.get(one_letter, 'ALA')

def get_active_site_indices_bio(pdb_file, cutoff=5.0):
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("temp", pdb_file)
    model = structure[0]
    
    metal_atoms = []
    target_resnames = ["MN", "FE", "HEM", "HEME"]
    
    for chain in model:
        for residue in chain:
            if residue.resname.strip().upper() in target_resnames:
                for atom in residue:
                    metal_atoms.append(atom)
    
    if not metal_atoms:
        return []

    active_site_indices = set()
    
    for chain in model:
        for residue in chain:
            if residue.resname.strip().upper() in target_resnames + ["HOH", "WAT"]: 
                continue
            
            for atom in residue:
                for metal in metal_atoms:
                    diff = atom.coord - metal.coord
                    dist = np.sqrt(np.sum(diff * diff))
                    if dist < cutoff:
                        active_site_indices.add(residue.id[1])
                        break
                if residue.id[1] in active_site_indices: break
                
    return list(active_site_indices)

def calculate_active_site_rmsd(pose1, pose2, pdb_ids):
    pose_res_list = []
    for pdb_id in pdb_ids:
        r = pose1.pdb_info().pdb2pose("A", pdb_id)
        if r != 0: 
            pose_res_list.append(r)
    
    if not pose_res_list: 
        return 99.9

    atom_map = map_core_id_AtomID_core_id_AtomID()
    for res in pose_res_list:
        atom_name = "CA"
        if not pose1.residue(res).has(atom_name) or not pose2.residue(res).has(atom_name):
            continue
            
        id1 = AtomID(pose1.residue(res).atom_index(atom_name), res)
        id2 = AtomID(pose2.residue(res).atom_index(atom_name), res)
        atom_map[id1] = id2
    
    res_vector = vector1_unsigned_long()
    for res in pose_res_list:
        res_vector.append(res)
    
    try:
        return rms_at_corresponding_atoms(pose1, pose2, atom_map, res_vector)
    except:
        return 99.9

def run_gen_task(args):
    combo_name, mutations, replicate_id = args
    init_rosetta()

    try:
        pose = pose_from_pdb(PDB_FILE)
        sfxn = create_score_function("ref2015_cart")

        relax = FastRelax()
        relax.set_scorefxn(sfxn)
        relax.cartesian(True)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
        relax.min_type("lbfgs_armijo_nonmonotone")

        wt_pose = pose.clone()
        relax.apply(wt_pose)
        wt_score = sfxn(wt_pose)
        
        wt_pdb_name = f"{combo_name}_WT_rep{replicate_id}.pdb"
        wt_pose.dump_pdb(wt_pdb_name)

        mut_pose = wt_pose.clone()
        for (old_res, res_num, new_res) in mutations:
            pose_res_id = pose.pdb_info().pdb2pose("A", res_num)
            
            if pose_res_id == 0:
                print(f"WARNING: Residue {res_num} not found in PDB!")
                continue
                
            target_res_3 = get_3_letter_code(new_res)
            MutateResidue(pose_res_id, target_res_3).apply(mut_pose)
        
        relax.apply(mut_pose)
        mut_score = sfxn(mut_pose)

        ddg = mut_score - wt_score
        
        active_indices = get_active_site_indices_bio(wt_pdb_name, ACTIVE_SITE_CUTOFF)
        rmsd = calculate_active_site_rmsd(wt_pose, mut_pose, active_indices)
        
        mut_pdb_name = f"{combo_name}_rep{replicate_id}.pdb"
        mut_pose.dump_pdb(mut_pdb_name)
        
        classification = "NON-FUNCTIONAL"
        if rmsd <= 0.5:
            classification = "OPTIMAL"
        elif rmsd <= 2.0:
            classification = "PERMISSIBLE"

        return (combo_name, replicate_id, ddg, rmsd, classification, mut_pdb_name)

    except Exception as e:
        print(f"Error in {combo_name} rep {replicate_id}: {e}")
        return None

def main():
    if not os.path.exists(PDB_FILE):
        print(f"ERROR: {PDB_FILE} not found")
        return
    
    multiprocessing.set_start_method("spawn", force=True)
    
    job_queue = []
    for name, muts in COMBOS.items():
        for i in range(1, N_REPLICATES + 1):
            job_queue.append((name, muts, i))

    print(f"Starting Final PDB Generation - {len(job_queue)} Total Jobs")
    print("-" * 80)
    print(f"{'Variant':<15} | {'Rep':<3} | {'ddG':<8} | {'RMSD':<8} | {'Class'}")

    results = []
    with multiprocessing.Pool(processes=THREADS_TO_USE) as pool:
        for result in pool.imap_unordered(run_gen_task, job_queue):
            if result:
                name, rep, ddg, rmsd, cls, fname = result
                print(f"{name:<15} | {rep:<3} | {ddg:<8.3f} | {rmsd:<8.3f} | {cls}")
                results.append(result)

    if not results:
        print("ERROR: No results generated")
        return

    df = pd.DataFrame(results, columns=["Variant", "Replicate", "ddG", "Active_Site_RMSD", "Classification", "PDB_File"])
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nDone! PDBs generated and scores saved to {OUTPUT_CSV}")

    print("\n" + "="*80)
    print("SELECTING REPRESENTATIVE PDBs (Closest to IQM)")
    print("="*80)

    for variant in COMBOS.keys():
        var_df = df[df['Variant'] == variant]
        if var_df.
