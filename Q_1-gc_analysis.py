"""
gc_analysis.py — Part A Q1: GC content analysis (fresh code).

Implements:
- overall_gc(seq)
- gc_sliding(seq, window, step)
- gc_skew_array(seq) : (G-C)/(G+C) per position
- compare_gc(seqs: dict[name]->seq)
- CLI with simple CSV export (no external libs required)

Run:
  python gc_analysis.py
"""

from __future__ import annotations
from typing import Dict, List, Tuple
from core_dna import clean_dna, validate_dna, sliding_windows

def _count_gc(s: str) -> int:
    return s.count("G") + s.count("C")

def overall_gc(seq: str) -> float:
    """Return GC% (0..100)."""
    seq = clean_dna(seq)
    n = len(seq)
    if n == 0: return 0.0
    return _count_gc(seq) * 100.0 / n

def gc_sliding(seq: str, window: int, step: int = 1) -> List[Tuple[int, float]]:
    """Return list of (start_index, GC%) over sliding windows."""
    seq = clean_dna(seq)
    out: List[Tuple[int, float]] = []
    for i, w in sliding_windows(seq, window, step):
        if len(w) < window: break
        out.append((i, _count_gc(w) * 100.0 / len(w)))
    return out

def gc_skew_array(seq: str) -> List[float]:
    """
    GC skew at each position i (cumulative):
    skew = (G - C) / (G + C) over prefix (avoid div by zero -> 0.0)
    """
    seq = clean_dna(seq)
    g = c = 0
    out: List[float] = []
    for ch in seq:
        if ch == "G": g += 1
        elif ch == "C": c += 1
        denom = g + c
        out.append(((g - c) / denom) if denom else 0.0)
    return out

def compare_gc(seqs: Dict[str, str]) -> Dict[str, float]:
    """Return {name: GC%} for multiple sequences."""
    res = {}
    for name, s in seqs.items():
        res[name] = overall_gc(s)
    return res

def _save_csv(rows: List[Tuple], path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")

# ---------------- CLI ----------------
if __name__ == "__main__":
    print("GC Analysis (Q1). Examples:")
    print("  overall GC, sliding-window GC, GC skew, compare multiple\n")
    while True:
        print("Menu:")
        print("  1) Overall GC% of a sequence")
        print("  2) Sliding-window GC%")
        print("  3) GC skew (cumulative) -> CSV")
        print("  4) Compare GC% across multiple sequences")
        print("  0) Exit")
        ch = input("Enter 1/2/3/4/0: ").strip()

        if ch == "0":
            break
        elif ch == "1":
            s = input("Enter sequence: ")
            cs = clean_dna(s)
            print(f"GC%: {overall_gc(cs):.2f}\n")
        elif ch == "2":
            s = input("Enter sequence: ")
            w = int(input("Window size (e.g., 100): ") or "100")
            step = int(input("Step (e.g., 10): ") or "10")
            cs = clean_dna(s)
            data = gc_sliding(cs, w, step)
            print(f"Computed {len(data)} windows. First 5:", data[:5], "\n")
            if input("Save CSV? (y/n): ").lower().startswith("y"):
                _save_csv(data, "gc_sliding.csv")
                print("Saved gc_sliding.csv\n")
        elif ch == "3":
            s = input("Enter sequence: ")
            cs = clean_dna(s)
            arr = gc_skew_array(cs)
            rows = [(i, val) for i, val in enumerate(arr)]
            _save_csv(rows, "gc_skew.csv")
            print("Saved gc_skew.csv (columns: index, skew)\n")
        elif ch == "4":
            n = int(input("How many sequences?: ") or "2")
            seqs = {}
            for i in range(1, n+1):
                name = input(f"Name {i}: ").strip() or f"S{i}"
                s = clean_dna(input(f"Sequence {i}: "))
                seqs[name] = s
            res = compare_gc(seqs)
            for k,v in res.items():
                print(f"{k}: {v:.2f}%")
            print()
        else:
            print("Enter valid option.\n")
