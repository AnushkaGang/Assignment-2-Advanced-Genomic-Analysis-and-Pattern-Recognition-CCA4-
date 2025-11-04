"""
composition.py — Part A Q2: Sequence composition analysis (fresh code).

Implements:
- di_tri_frequencies(seq)
- find_cpg_islands(seq, window=200, min_gc=50.0, min_oe=0.6)
- codon_usage(seq)  (expects coding DNA, frame 0)
- classify_by_composition(seq) -> simple label

Run:
  python composition.py
"""

from __future__ import annotations
from typing import Dict, List, Tuple
from core_dna import clean_dna, sliding_windows

def di_tri_frequencies(seq: str) -> Dict[str, float]:
    """Return normalized di- and tri-nucleotide frequencies (0..1)."""
    s = clean_dna(seq)
    n = len(s)
    res: Dict[str, int] = {}
    # dinucleotides
    for i in range(n-1):
        k = s[i:i+2]
        res[k] = res.get(k, 0) + 1
    # trinucleotides
    for i in range(n-2):
        k = s[i:i+3]
        res[k] = res.get(k, 0) + 1
    total = sum(res.values()) or 1
    return {k: v/total for k, v in res.items()}

def _obs_exp_cpg(window: str) -> float:
    """Observed/Expected CpG ratio in a window."""
    g = window.count("G")
    c = window.count("C")
    cg = sum(1 for i in range(len(window)-1) if window[i:i+2] == "CG")
    exp = (c * g) / max(len(window),1)
    return (cg / exp) if exp > 0 else 0.0

def find_cpg_islands(seq: str, window: int = 200,
                     min_gc: float = 50.0, min_oe: float = 0.6) -> List[Tuple[int,int,float,float]]:
    """
    Identify CpG-like regions using common heuristic:
      - GC% >= min_gc
      - Observed/Expected(CpG) >= min_oe
    Returns list of (start, end, gc%, O/E).
    """
    s = clean_dna(seq)
    out: List[Tuple[int,int,float,float]] = []
    for i, w in sliding_windows(s, window, step=1):
        gc = (w.count("G")+w.count("C"))*100.0/len(w)
        oe = _obs_exp_cpg(w)
        if gc >= min_gc and oe >= min_oe:
            out.append((i, i+window, gc, oe))
    return out

_CODON_TABLE = {
    # U is not used; we assume DNA (T)
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def codon_usage(seq: str) -> Dict[str, float]:
    """Return frequency distribution of codons in frame 0 (ignore last <3)."""
    s = clean_dna(seq)
    counts: Dict[str,int] = {}
    for i in range(0, len(s) - 2, 3):
        cod = s[i:i+3]
        if len(cod) == 3:
            counts[cod] = counts.get(cod,0)+1
    total = sum(counts.values()) or 1
    return {k: v/total for k,v in counts.items()}

def classify_by_composition(seq: str) -> str:
    """Very simple classifier demo: GC-rich/AT-rich/neutral."""
    s = clean_dna(seq)
    n = len(s) or 1
    gc = (s.count("G")+s.count("C"))*100.0/n
    if gc >= 60: return "GC-rich"
    if gc <= 40: return "AT-rich"
    return "balanced"

# --------------- CLI ---------------
if __name__ == "__main__":
    print("Composition (Q2).")
    while True:
        print(" 1) di/tri-nucleotide frequencies")
        print(" 2) find CpG islands")
        print(" 3) codon usage (frame 0)")
        print(" 4) simple composition classification")
        print(" 0) Exit")
        ch = input("Enter: ").strip()
        if ch == "0": break
        if ch == "1":
            s = input("Sequence: ")
            print(di_tri_frequencies(s)); print()
        elif ch == "2":
            s = input("Sequence: ")
            w = int(input("window [200]: ") or "200")
            gc = float(input("min GC% [50]: ") or "50")
            oe = float(input("min O/E [0.6]: ") or "0.6")
            out = find_cpg_islands(s, w, gc, oe)
            print("start,end,GC%,O/E:")
            for r in out[:20]:
                print(r)
            print(f"Total islands: {len(out)}\n")
        elif ch == "3":
            s = input("Coding DNA (frame 0): ")
            print(codon_usage(s)); print()
        elif ch == "4":
            s = input("Sequence: ")
            print("Class:", classify_by_composition(s), "\n")
        else:
            print("valid options are 0..4\n")
