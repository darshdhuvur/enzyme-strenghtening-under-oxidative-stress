]import os
import glob
import numpy as np
from collections import Counter, defaultdict

OUTPUT_FOLDER = "mpnn_1DGB_output/seqs"
SCORE_FOLDER = "mpnn_1DGB_output/scores"

def parse_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    entries = []
    current_header = ""
    current_seq = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                entries.append((current_header, current_seq))
            current_header = line
            current_seq = ""
        else:
            current_seq += line
    if current_header:
        entries.append((current_header, current_seq))
        
    return entries

def get_mutations(wt_seq, des_seq):
    muts = []
    for i, (w, d) in enumerate(zip(wt_seq, des_seq)):
        if w != d:
            muts.append(f"{w}{i+1}{d}")
    return muts

def load_scores():
    score_files = glob.glob(os.path.join(SCORE_FOLDER, "*.npz"))
    if not score_files:
        return None
    
    try:
        scores_data = np.load(score_files[0])
        scores = scores_data['score']
        print(f"Loaded {len(scores)} scores from {os.path.basename(score_files[0])}")
        print(f"Score range: {scores.min():.3f} to {scores.max():.3f}")
        print(f"Mean score: {scores.mean():.3f}")
        return scores
    except Exception as e:
        print(f"Warning: Could not load scores: {e}")
        return None

def analyze_mutations(sequences, wt_seq, scores=None, top_percent=None):
    if scores is not None and top_percent is not None:
        n_top = int(len(sequences) * top_percent / 100)
        sorted_indices = np.argsort(scores)[:n_top]
        sequences = [sequences[i] for i in sorted_indices]
        print(f"\nFiltered to top {top_percent}% by score ({len(sequences)} sequences)")
        print(f"Best score: {scores[sorted_indices[0]]:.3f}")
        print(f"Worst score in selection: {scores[sorted_indices[-1]]:.3f}")
    
    mutation_counter = Counter()
    mutation_scores = defaultdict(list)
    
    for idx, (header, seq) in enumerate(sequences):
        muts = get_mutations(wt_seq, seq)
        mutation_counter.update(muts)
        
        if scores is not None:
            for mut in muts:
                mutation_scores[mut].append(scores[idx])
    
    return mutation_counter, mutation_scores, len(sequences)

def rank_mutations(mutation_counter, mutation_scores, total_sequences, min_frequency=5):
    ranked = []
    
    for mut, count in mutation_counter.items():
        if count < min_frequency:
            continue
            
        frequency_pct = (count / total_sequences) * 100
        
        avg_score = np.mean(mutation_scores[mut]) if mut in mutation_scores else 0
        
        if mutation_scores:
            rank_score = avg_score - (frequency_pct / 100)
        else:
            rank_score = -frequency_pct
        
        ranked.append((mut, count, frequency_pct, avg_score, rank_score))
    
    ranked.sort(key=lambda x: x[4])
    
    return ranked

def main():
    if not os.path.exists(OUTPUT_FOLDER):
        print(f"Error: Folder '{OUTPUT_FOLDER}' not found.")
        print("Did you run the ProteinMPNN command successfully?")
        return

    fa_files = glob.glob(os.path.join(OUTPUT_FOLDER, "*.fa"))
    if not fa_files:
        print(f"No .fa files found in {OUTPUT_FOLDER}")
        return

    print("="*80)
    print("PROTEINMPNN MUTATION ANALYSIS")
    print("="*80)
    
    scores = load_scores()
    
    first_entries = parse_fasta(fa_files[0])
    wt_header, wt_seq = first_entries[0]
    print(f"\nWildtype sequence length: {len(wt_seq)}")
    
    all_sequences = []
    for fa_file in fa_files:
        entries = parse_fasta(fa_file)
        designs = entries[1:]
        all_sequences.extend(designs)
    
    print(f"Total designed sequences: {len(all_sequences)}")
    
    print("\n" + "="*80)
    print("ANALYSIS 1: ALL SEQUENCES")
    print("="*80)
    
    mut_counter_all, mut_scores_all, total_all = analyze_mutations(
        all_sequences, wt_seq, scores=scores
    )
    ranked_all = rank_mutations(mut_counter_all, mut_scores_all, total_all, min_frequency=5)
    
    print(f"\nTop 20 mutations (all sequences)")
    print(f"{'Rank':<6}{'Mutation':<12}{'Count':<8}{'Freq %':<10}{'Avg Score':<12}{'Rank Score':<12}")
    print("-" * 80)
    for i, (mut, count, freq, avg_score, rank_score) in enumerate(ranked_all[:20], 1):
        score_str = f"{avg_score:.3f}" if scores is not None else "N/A"
        rank_str = f"{rank_score:.3f}" if scores is not None else f"{rank_score:.3f}"
        print(f"{i:<6}{mut:<12}{count:<8}{freq:<10.1f}{score_str:<12}{rank_str:<12}")
    
    if scores is not None:
        print("\n" + "="*80)
        print("ANALYSIS 2: TOP 20% BY SCORE (HIGHEST CONFIDENCE)")
        print("="*80)
        
        mut_counter_top, mut_scores_top, total_top = analyze_mutations(
            all_sequences, wt_seq, scores=scores, top_percent=20
        )
        ranked_top = rank_mutations(mut_counter_top, mut_scores_top, total_top, min_frequency=2)
        
        print(f"\nTop 20 mutations (in best 20% of designs)")
        print(f"{'Rank':<6}{'Mutation':<12}{'Count':<8}{'Freq %':<10}{'Avg Score':<12}{'Rank Score':<12}")
        print("-" * 80)
        for i, (mut, count, freq, avg_score, rank_score) in enumerate(ranked_top[:20], 1):
            print(f"{i:<6}{mut:<12}{count:<8}{freq:<10.1f}{avg_score:<12.3f}{rank_score:<12.3f}")
    
    print("\n" + "="*80)
    print("RECOMMENDATIONS FOR ROSETTA FASTRELAX")
    print("="*80)
    
    if scores is not None:
        print("\nUsing score-filtered mutations (recommended)")
        print(f"\nSuggested mutations to test (appear in >30% of top 20% designs):")
        recommendations = [
            (mut, count, freq, avg_score) 
            for mut, count, freq, avg_score, _ in ranked_top 
            if freq > 30
        ]
    else:
        print("\nNo scores available - using frequency only")
        print(f"\nSuggested mutations to test (appear in >50% of all designs):")
        recommendations = [
            (mut, count, freq, avg_score) 
            for mut, count, freq, avg_score, _ in ranked_all 
            if freq > 50
        ]
    
    if recommendations:
        print(f"\n{'Mutation':<12}{'Count':<8}{'Frequency':<12}{'Avg Score':<12}")
        print("-" * 50)
        for mut, count, freq, avg_score in recommendations[:10]:
            score_str = f"{avg_score:.3f}" if scores is not None else "N/A"
            print(f"{mut:<12}{count:<8}{freq:<12.1f}%{score_str:<12}")
        
        print(f"\nRecommended: Test these {len(recommendations[:10])} mutations + 1 WT control")
        print(f"Total FastRelax runs: {len(recommendations[:10]) + 1} variants × 25 runs = {(len(recommendations[:10]) + 1) * 25} jobs")
    else:
        print("\nNo mutations meet the frequency threshold")
        print("Consider lowering threshold or generating more sequences")
    
    output_file = "mutation_1DGB_analysis_results.txt"
    with open(output_file, 'w') as f:
        f.write("COMPLETE MUTATION RANKING\n")
        f.write("="*80 + "\n\n")
        f.write(f"{'Rank':<6}{'Mutation':<12}{'Count':<8}{'Freq %':<10}{'Avg Score':<12}\n")
        f.write("-" * 80 + "\n")
        
        results_to_save = ranked_top if scores is not None else ranked_all
        for i, (mut, count, freq, avg_score, _) in enumerate(results_to_save, 1):
            score_str = f"{avg_score:.3f}" if scores is not None else "N/A"
            f.write(f"{i:<6}{mut:<12}{count:<8}{freq:<10.1f}{score_str:<12}\n")
    
    print(f"\nFull results saved to: {output_file}")

if __name__ == "__main__":
    main()
