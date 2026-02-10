import pyrosetta
from pyrosetta import init, pose_from_pdb, create_score_function
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
import multiprocessing
import os
import sys

PDB_FILE = "5T30_cleaned_monomer.pdb"
RAW_OUTPUT_FILE = "cartesian_ddg_t10_5T30_raw.csv"
THREADS_TO_USE = 14
N_REPLICATES = 25

mutations_to_run = [
    "D10P",
    "G12D",
    "K29V",
    "S75T",
    "K106A",
    "S115A",
    "N129D",
    "C140T",
    "K170G",
    "V172D"
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

def calculate_ddg_task(args):
    mutation_str, replicate_id = args
    init_rosetta()

    try:
        old_res = mutation_str[0]
        new_res = mutation_str[-1]
        res_num = int(mutation_str[1:-1])

        pose = pose_from_pdb(PDB_FILE)
        
        pose_res_id = pose.pdb_info().pdb2pose("A", res_num)

        if pose_res_id == 0:
            return (mutation_str, replicate_id, 999999.0, 999999.0, 999999.0)

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

        ddg_value = mut_score - wt_score

        return (mutation_str, replicate_id, ddg_value, wt_score, mut_score)

    except Exception:
        return (mutation_str, replicate_id, 999999.0, 999999.0, 999999.0)

def main():
    multiprocessing.set_start_method("spawn", force=True)

    completed_jobs = set()
    if os.path.exists(RAW_OUTPUT_FILE):
        with open(RAW_OUTPUT_FILE, "r") as f:
            for line in f:
                if line.startswith("Mutation"):
                    continue
                parts = line.strip().split(",")
                if len(parts) >= 2:
                    completed_jobs.add((parts[0], int(parts[1])))
    else:
        with open(RAW_OUTPUT_FILE, "w") as f:
            f.write("Mutation,Replicate,ddG,WT_Score,Mut_Score\n")

    job_queue = []
    for mut in mutations_to_run:
        for i in range(N_REPLICATES):
            replicate_id = i + 1
            if (mut, replicate_id) not in completed_jobs:
                job_queue.append((mut, replicate_id))

    if not job_queue:
        print("All jobs completed.")
    else:
        print(f"Starting processing of {len(job_queue)} jobs...")

        with multiprocessing.Pool(processes=THREADS_TO_USE) as pool:
            for result in pool.imap_unordered(calculate_ddg_task, job_queue):
                mut, rep, ddg, wt, mt = result
                
                print(f"Finished: {mut} (Rep {rep}) -> ddG: {ddg:.2f}")

                with open(RAW_OUTPUT_FILE, "a") as f:
                    f.write(f"{mut},{rep},{ddg},{wt},{mt}\n")

if __name__ == "__main__":
    main()
