import pyrosetta
from pyrosetta import init, pose_from_pdb, create_score_function
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.core.scoring import rms_at_corresponding_atoms, CA_rmsd
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID
from pyrosetta.rosetta.utility import vector1_unsigned_long
import multiprocessing
import os
import pandas as pd
import shutil
import numpy as np

PDB_FILE = "5T30_cleaned_monomer.pdb"
OUTPUT_CSV = "methionine_replacement_results_5T30.csv"
THREADS_TO_USE = 15
N_REPLICATES = 25
ACTIVE_SITE_PDB_IDS = [26, 74, 159, 163]

MUTATIONS_TO_RUN = [
    "M23L",
    "M194L"
]

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
        try:
            return rms_at_corresponding_atoms(pose1, pose2, atom_map)
        except:
            return CA_rmsd(pose1, pose2)

def validate_task(args):
    mutation_str, replicate_id = args
    init_rosetta()

    try:
        old_res = mutation_str[0]
        new_res = mutation_str[-1]
        res_num = int(mutation_str[1:-1])

        pose = pose_from_pdb(PDB_FILE)
        pose_res_id = pose.pdb_info().pdb2pose("A", res_num)

        if pose_res_id == 0:
            return None

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

        mut_pose = pose.clone()
        target_res_3 = get_3_letter_code(new_res)
        MutateResidue(pose_res_id, target_res_3).apply(mut_pose)
        relax.apply(mut_pose)
        mut_score = sfxn(mut_pose)

        ddg = mut_score - wt_score
        rmsd = calculate_active_site_rmsd(wt_pose, mut_pose, ACTIVE_SITE_PDB_IDS)
        full_rmsd = CA_rmsd(wt_pose, mut_pose)

        pdb_filename = f"{mutation_str}_rep{replicate_id}.pdb"
        mut_pose.dump_pdb(pdb_filename)

        return (mutation_str, replicate_id, ddg, rmsd, full_rmsd, pdb_filename)

    except Exception as e:
        print(f"Error in {mutation_str} rep {replicate_id}: {e}")
        return None

def calculate_iqm(values):
    if not values: return 0.0
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    filtered = [v for v in values if q1 <= v <= q3]
    if not filtered: return np.mean(values)
    return np.mean(filtered)

def main():
    if not os.path.exists(PDB_FILE):
        print(f"ERROR: {PDB_FILE} not found")
        return
    
    multiprocessing.set_start_method("spawn", force=True)
    
    job_queue = []
    for mut in MUTATIONS_TO_RUN:
        for i in range(N_REPLICATES):
            job_queue.append((mut, i + 1))

    print(f"Starting validation of {len(job_queue)} jobs...")
    print(f"{'Mutation':<10} | {'Rep':<3} | {'ddG':<8} | {'ActiveRMSD':<10} | {'FullRMSD':<8} | {'File'}")
    print("-" * 80)

    results = []
    with multiprocessing.Pool(processes=THREADS_TO_USE) as pool:
        for result in pool.imap_unordered(validate_task, job_queue):
            if result:
                mut, rep, ddg, rmsd, full_rmsd, fname = result
                print(f"{mut:<10} | {rep:<3} | {ddg:<8.3f} | {rmsd:<10.3f} | {full_rmsd:<8.3f} | {fname}")
                results.append(result)

    if not results:
        print("ERROR: No results generated")
        return

    df = pd.DataFrame(results, columns=["Mutation", "Replicate", "ddG", "Active_Site_RMSD", "Full_RMSD", "PDB_File"])
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nResults saved to {OUTPUT_CSV}")
    
    print("\n" + "="*80)
    print("SELECTING REPRESENTATIVE PDBs (Closest to IQM ddG)")
    print("="*80)
    
    for mut in MUTATIONS_TO_RUN:
        mut_df = df[df['Mutation'] == mut].copy()
        if mut_df.empty:
            print(f"WARNING: No results for {mut}")
            continue
            
        ddg_values = mut_df['ddG'].tolist()
        iqm_score = calculate_iqm(ddg_values)
        
        mut_df['Diff_to_IQM'] = abs(mut_df['ddG'] - iqm_score)
        best_row = mut_df.loc[mut_df['Diff_to_IQM'].idxmin()]
        
        best_file = best_row['PDB_File']
        best_ddg = best_row['ddG']
        best_rmsd = best_row['Active_Site_RMSD']
        best_full_rmsd = best_row['Full_RMSD']
        
        final_name = f"{mut}_Representative.pdb"
        
        if os.path.exists(best_file):
            shutil.copy(best_file, final_name)
            
            print(f"\n{mut}:")
            print(f"  Calculated IQM ddG: {iqm_score:7.3f}")
            print(f"  Selected Rep:       {best_row['Replicate']} (ddG: {best_ddg:7.3f})")
            print(f"  Active Site RMSD:   {best_rmsd:7.3f} A")
            print(f"  Full RMSD:          {best_full_rmsd:7.3f} A")
            print(f"  Saved As:           {final_name}")
            
            if best_rmsd > 1.0:
                print(f"  WARNING: High active site RMSD - possible active site distortion")
            elif best_rmsd > 0.5:
                print(f"  CAUTION: Moderate active site RMSD - check structure in YASARA")
            else:
                print(f"  OK: Low active site RMSD - likely safe")
        else:
            print(f"WARNING: {best_file} not found")

        for f in mut_df['PDB_File']:
            if f != best_file and os.path.exists(f):
                os.remove(f)
        if os.path.exists(best_file) and best_file != final_name:
             os.remove(best_file)
    
    print("\n" + "="*80)
    print(f"Total structures processed: {len(results)}")
    print(f"Representative structures: {len([f for f in os.listdir('.') if f.endswith('_Representative.pdb')])}")

if __name__ == "__main__":
    main()
