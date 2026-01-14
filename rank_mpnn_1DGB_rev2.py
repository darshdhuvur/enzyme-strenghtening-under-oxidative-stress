import os
from collections import Counter

output_folder = "mpnn_jan13_1DGB_output/seqs"

if not os.path.exists(output_folder):
    print(f"Error: Folder {output_folder} not found.")
    exit()

fa_files = [f for f in os.listdir(output_folder) if f.endswith('.fa')]
if not fa_files:
    print(f"No .fa files found in {output_folder}")
    exit()

target_file = os.path.join(output_folder, fa_files[0])
print(f"Analyzing {target_file}...")

def get_mutations(wt_seq, des_seq):
    muts = []
    for i, (w, d) in enumerate(zip(wt_seq, des_seq)):
        if w != d:
            muts.append(f"{w}{i+1}{d}")
    return muts

with open(target_file, 'r') as f:
    lines = f.readlines()

sequences = []
current_seq = ""
for line in lines:
    line = line.strip()
    if line.startswith(">"):
        if current_seq:
            sequences.append(current_seq)
            current_seq = ""
    else:
        current_seq += line
if current_seq:
    sequences.append(current_seq)

wt_seq = sequences[0]
designs = sequences[1:]

all_mutations = []
for des in designs:
    muts = get_mutations(wt_seq, des)
    all_mutations.extend(muts)

counts = Counter(all_mutations)
print("\nTop 10 Most Frequent Mutations:")
print("-" * 30)
for mut, count in counts.most_common(10):
    print(f"{mut}: {count}/{len(designs)} sequences")
print("-" * 30)
