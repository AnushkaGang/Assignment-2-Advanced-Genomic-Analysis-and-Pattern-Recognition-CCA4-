"""
test_all.py — CCA4 Assignment Testing
This file has basic functional tests for each module.
Run:  python test_all.py
"""

from gc_analysis import overall_gc
from composition import classify_by_composition
from hamming_distance import hamming_distance
from motif_finding import kmp_search, approx_matches
from orf_reading_frames import find_all_orfs
from protein_translation import translate_dna

def main():
    print("\n--- BASIC TEST SUITE ---\n")

    seq = "ATGCGCGATATATATGCGC"

    print("GC%:", overall_gc(seq))

    print("Composition label:", classify_by_composition(seq))

    print("Hamming test (same):", hamming_distance("AAAA","AAAA"))
    print("Hamming test (diff):", hamming_distance("AAAA","AAAT"))

    print("KMP search:", kmp_search(seq, "GCG"))
    print("Approx search (1 mismatch):", approx_matches(seq, "GCG", 1))

    orfs = find_all_orfs(seq, min_len_codon=2)
    print("ORFs found:", sum(len(v) for v in orfs.values()))

    print("Translate F0:", translate_dna(seq, frame=0)[:20])

    print("\n--- TESTS DONE ---\n")

if __name__ == "__main__":
    main()
