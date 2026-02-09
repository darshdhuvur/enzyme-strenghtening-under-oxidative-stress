import pyrosetta
from pyrosetta import init, pose_from_pdb, create_score_function
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
import multiprocessing
import os
import pandas as pd

PDB_FILE = "1DGB_cleaned_monomer.pdb"
OUTPUT_CSV = "final_1DGB_mutants_results.csv"
THREADS_TO_USE = 12
N_REPLICATES = 25

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

        mut_pose = pose.clone()
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

        fname = f"{combo_name}_rep{replicate_id}.pdb"
        mut_pose.dump_pdb(fname)

        return (combo_name, replicate_id, ddg, fname)

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
    print("-" * 60)
    print(f"{'Variant':<15} | {'Rep':<3} | {'ddG':<8} | {'File'}")

    results = []
    with multiprocessing.Pool(processes=THREADS_TO_USE) as pool:
        for result in pool.imap_unordered(run_gen_task, job_queue):
            if result:
                name, rep, ddg, fname = result
                print(f"{name:<15} | {rep:<3} | {ddg:<8.3f} | {fname}")
                results.append(result)

    if not results:
        print("ERROR: No results generated")
        return

    df = pd.DataFrame(results, columns=["Variant", "Replicate", "ddG", "PDB_File"])
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nDone! PDBs generated and scores saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
