"""
motif_finding.py — Part B Q4: Motif finding & pattern matching (fresh code).

Implements:
- kmp_search(text, pattern) -> all start indices (exact)
- approx_matches(text, pattern, k) -> allow up to k mismatches (Hamming)
- find_all_occurrences(seqs: dict[name]->text, pattern, mismatches)
- conservation_score(seqs, pattern, mismatches) -> fraction sequences containing it

Run:
  python motif_finding.py
"""

from __future__ import annotations
from typing import List, Dict, Tuple
from core_dna import clean_dna

def _kmp_prefix(p: str) -> List[int]:
    pi = [0]*len(p)
    j = 0
    for i in range(1, len(p)):
        while j > 0 and p[i] != p[j]:
            j = pi[j-1]
        if p[i] == p[j]:
            j += 1
            pi[i] = j
    return pi

def kmp_search(text: str, pattern: str) -> List[int]:
    """Return all indices where pattern occurs in text (exact)."""
    t = clean_dna(text); p = clean_dna(pattern)
    if not p: return []
    pi = _kmp_prefix(p)
    out: List[int] = []
    j = 0
    for i in range(len(t)):
        while j > 0 and t[i] != p[j]:
            j = pi[j-1]
        if t[i] == p[j]:
            j += 1
            if j == len(p):
                out.append(i - j + 1)
                j = pi[j-1]
    return out

def approx_matches(text: str, pattern: str, k: int) -> List[int]:
    """Return all start indices with Hamming distance <= k."""
    t = clean_dna(text); p = clean_dna(pattern)
    m, n = len(p), len(t)
    if m == 0 or m > n: return []
    out: List[int] = []
    for i in range(n - m + 1):
        window = t[i:i+m]
        d = sum(1 for a,b in zip(window, p) if a != b)
        if d <= k:
            out.append(i)
    return out

def find_all_occurrences(seqs: Dict[str, str], pattern: str, mismatches: int = 0) -> Dict[str, List[int]]:
    res: Dict[str, List[int]] = {}
    for name, s in seqs.items():
        if mismatches == 0:
            res[name] = kmp_search(s, pattern)
        else:
            res[name] = approx_matches(s, pattern, mismatches)
    return res

def conservation_score(seqs: Dict[str, str], pattern: str, mismatches: int = 0) -> float:
    """Fraction of sequences where pattern occurs at least once (with <= mismatches)."""
    occ = find_all_occurrences(seqs, pattern, mismatches)
    hits = sum(1 for v in occ.values() if len(v) > 0)
    total = len(seqs) or 1
    return hits / total

# --------------- CLI ---------------
if __name__ == "__main__":
    print("Motif Finding (Q4)")
    t = input("Enter main sequence (or leave blank for multi-seq mode): ").strip()
    pat = input("Pattern: ").strip()
    k = int(input("Allowed mismatches [0]: ") or "0")
    if t:
        inds = approx_matches(t, pat, k) if k>0 else kmp_search(t, pat)
        print("Indices:", inds)
    else:
        n = int(input("How many sequences [3]: ") or "3")
        seqs = {}
        for i in range(1, n+1):
            name = input(f"Name {i}: ").strip() or f"S{i}"
            seqs[name] = input(f"Seq {i}: ")
        occ = find_all_occurrences(seqs, pat, k)
        for name, idxs in occ.items():
            print(f"{name}: {idxs}")
        print("Conservation score:", conservation_score(seqs, pat, k))
