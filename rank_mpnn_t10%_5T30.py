import os
import glob
import numpy as np
from collections import Counter

OUTPUT_FOLDER = "mpnn_5T30_output"
SEQ_FOLDER = os.path.join(OUTPUT_FOLDER, "seqs")
TOP_PERCENT = 0.10

def parse_fasta(file_path):
    entries = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    current_header = None
    current_seq = []
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                entries.append({"header": current_header, "seq": "".join(current_seq)})
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_header:
        entries.append({"header": current_header, "seq": "".join(current_seq)})
        
    return entries

def get_score(header):
    try:
        parts = header.split(',')
        for part in parts:
            if "score=" in part:
                return float(part.split('=')[1])
    except:
        return 999.9
    return 999.9

def get_mutations(wt_seq, mut_seq):
    muts = []
    for i, (wt, mut) in enumerate(zip(wt_seq, mut_seq)):
        if wt != mut:
            muts.append(f"{wt}{i+1}{mut}")
    return muts

fasta_files = glob.glob(os.path.join(SEQ_FOLDER, "*.fa"))
if not fasta_files:
    print(f"Error: No .fa files found in {SEQ_FOLDER}")
    exit()

fasta_file = fasta_files[0]
print(f"Analyzing: {fasta_file}")

data = parse_fasta(fasta_file)

wt_entry = data[0]
wt_seq = wt_entry['seq']
designs = data[1:]

print(f"Wildtype Length: {len(wt_seq)}")
print(f"Total Designs: {len(designs)}")

scored_designs = []
for d in designs:
    score = get_score(d['header'])
    scored_designs.append({'seq': d['seq'], 'score': score})

scored_designs.sort(key=lambda x: x['score'])

cutoff_index = int(len(scored_designs) * TOP_PERCENT)
top_designs = scored_designs[:cutoff_index]

print(f"Filtering for Top {TOP_PERCENT*100}% Best Scoring Sequences")
print(f"Analyzing top {len(top_designs)} sequences...")
print(f"Best Score: {top_designs[0]['score']:.4f}")
print(f"Worst Score (in top 10%): {top_designs[-1]['score']:.4f}")

mutation_counts = Counter()
for d in top_designs:
    muts = get_mutations(wt_seq, d['seq'])
    mutation_counts.update(muts)

print("\n" + "="*60)
print(f"TOP MUTATIONS (In Top {TOP_PERCENT*100}% Best Scoring Designs)")
print("="*60)
print(f"{'Rank':<5} {'Mutation':<10} {'Count':<10} {'Freq %':<10}")
print("-" * 40)

sorted_muts = mutation_counts.most_common(30)
for rank, (mut, count) in enumerate(sorted_muts, 1):
    freq = (count / len(top_designs)) * 100
    print(f"{rank:<5} {mut:<10} {count:<10} {freq:.1f}")

print("\nDone.")
