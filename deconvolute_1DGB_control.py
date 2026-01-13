import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
import multiprocessing
import time
import os
import numpy as np
import pandas as pd

PDB_FILE = "1DGB_Relaxed_WT.pdb"
RAW_OUTPUT_FILE = "deconvolute_control_raw.csv"
FINAL_SUMMARY_FILE = "deconvolute_control_summary.csv"

THREADS_TO_USE = 6
N_STRUCT = 3

def init_rosetta():
    try:
        init("-ex1 -ex2 -use_input_sc -no_optH false -mute all")
    except:
        pass

def get_3_letter_code(one_letter):
    mapping = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    return mapping.get(one_letter, 'ALA')

def perform_mutation_task(args):
    mutation_str, run_id = args
    
    init_rosetta()
    
    old_res = mutation_str[0]
    new_res = mutation_str[-1]
    res_num = int(mutation_str[1:-1])
    
    target_res_3_letter = get_3_letter_code(new_res)
    
    pose = pose_from_pdb(PDB_FILE)
    pose_res_id = pose.pdb_info().pdb2pose("A", res_num)
    
    if pose_res_id == 0:
        return (mutation_str, run_id, 999999)

    pyrosetta.rosetta.protocols.simple_moves.MutateResidue(pose_res_id, target_res_3_letter).apply(pose)
    
    sfxn = create_score_function("ref2015_cart")
    relax = FastRelax()
    relax.set_scorefxn(sfxn)
    relax.cartesian(True)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)
    
    relax.apply(pose)
    
    final_score = sfxn(pose)
    
    return (mutation_str, run_id, final_score)

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    
    mutations = ["A5A"]
    
    completed_jobs = set()
    if os.path.exists(RAW_OUTPUT_FILE):
        print(f"Found existing data in {RAW_OUTPUT_FILE}. Checking for resume...")
        with open(RAW_OUTPUT_FILE, "r") as f:
            for line in f:
                if line.startswith("Mutation"): continue
                parts = line.strip().split(",")
                if len(parts) >= 2:
                    completed_jobs.add((parts[0], int(parts[1])))
    else:
        with open(RAW_OUTPUT_FILE, "w") as f:
            f.write("Mutation,RunID,Score\n")

    job_queue = []
    for mut in mutations:
        for i in range(N_STRUCT):
            run_id = i + 1
            if (mut, run_id) not in completed_jobs:
                job_queue.append((mut, run_id))
            else:
                print(f"Skipping {mut} Run {run_id}")
            
    print(f"Starting Control Run (A5A Baseline).")
    print(f"Total jobs remaining: {len(job_queue)}")
    
    start_time = time.time()
    
    if len(job_queue) > 0:
        with multiprocessing.Pool(processes=THREADS_TO_USE) as pool:
            for result in pool.imap_unordered(perform_mutation_task, job_queue):
                mut, run_id, score = result
                print(f"Finished {mut} Run {run_id}: {score:.2f}")
                
                with open(RAW_OUTPUT_FILE, "a") as f:
                    f.write(f"{mut},{run_id},{score}\n")
    
    print("-" * 40)
    print("Processing Final Stats...")
    
    if os.path.exists(RAW_OUTPUT_FILE):
        df = pd.read_csv(RAW_OUTPUT_FILE)
        summary = df.groupby("Mutation")["Score"].agg(["mean", "std", "count"])
        summary.to_csv(FINAL_SUMMARY_FILE)
        print(summary)
        print(f"\nSummary saved to {FINAL_SUMMARY_FILE}")
    
    print(f"Done in {time.time() - start_time:.2f} seconds.")
