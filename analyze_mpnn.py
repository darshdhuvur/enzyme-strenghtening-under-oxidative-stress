import os

fasta_file = "output_1GV3/seqs/1GV3_MASTER.fa"

if not os.path.exists(fasta_file):
    print(f"Error: Could not find {fasta_file}")
    exit()

print(f"Reading {fasta_file}...")

sequences = []
headers = []

with open(fasta_file, "r") as f:
    current_seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                sequences.append(current_seq)
            headers.append(line)
            current_seq = ""
        else:
            current_seq += line
    if current_seq:
        sequences.append(current_seq)

wt_seq = sequences[0]
print(f"Wild Type Length: {len(wt_seq)}")

mutations_found = set()

for i in range(1, len(sequences)):
    seq = sequences[i]
    for pos, (aa_wt, aa_new) in enumerate(zip(wt_seq, seq)):
        if aa_wt != aa_new:
            mutation = f"{aa_wt}{pos+1}{aa_new}"
            mutations_found.add(mutation)

print(f"\nFound {len(mutations_found)} unique suggested mutations.")
print("Top 10 suggestions:")
for m in list(mutations_found)[:10]:
    print(m)
