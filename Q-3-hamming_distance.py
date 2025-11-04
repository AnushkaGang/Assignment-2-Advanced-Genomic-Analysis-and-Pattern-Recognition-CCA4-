"""
hamming_distance.py — Part B Q3: Hamming distance & matrices (fresh code).

Implements:
- hamming_distance(a,b): counts mismatches + length difference
- distance_matrix(seqs: dict[name]->seq): symmetric matrix
- save_matrix_csv(mat, names, path)
- simple ASCII visualization (heatmap-like characters)

Run:
  python hamming_distance.py
"""

from __future__ import annotations
from typing import Dict, List, Tuple
from core_dna import clean_dna

def hamming_distance(a: str, b: str) -> int:
    """
    Hamming distance with unequal length handling:
      distance = mismatches_on_overlap + abs(len(a)-len(b))
    """
    a, b = clean_dna(a), clean_dna(b)
    n = min(len(a), len(b))
    d = sum(1 for i in range(n) if a[i] != b[i])
    d += abs(len(a) - len(b))
    return d

def distance_matrix(seqs: Dict[str, str]) -> Tuple[List[str], List[List[int]]]:
    names = list(seqs.keys())
    n = len(names)
    mat = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            d = hamming_distance(seqs[names[i]], seqs[names[j]])
            mat[i][j] = mat[j][i] = d
    return names, mat

def save_matrix_csv(names: List[str], mat: List[List[int]], path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("," + ",".join(names) + "\n")
        for name, row in zip(names, mat):
            f.write(name + "," + ",".join(map(str, row)) + "\n")

def ascii_visual(mat: List[List[int]]) -> str:
    """
    Simple ASCII visualization:
      small numbers -> '.'  medium -> '*'  large -> '#'
    """
    import math
    flat = [x for row in mat for x in row]
    mx = max(flat) if flat else 1
    lines = []
    for row in mat:
        line = ""
        for v in row:
            t = v / mx if mx else 0
            if t < 0.33: ch = "."
            elif t < 0.66: ch = "*"
            else: ch = "#"
            line += ch
        lines.append(line)
    return "\n".join(lines)

# --------------- CLI ---------------
if __name__ == "__main__":
    print("Hamming Distance (Q3)")
    n = int(input("How many sequences [3]: ") or "3")
    seqs = {}
    for i in range(1, n+1):
        name = input(f"Name {i}: ").strip() or f"S{i}"
        seqs[name] = input(f"Seq {i}: ")
    names, mat = distance_matrix(seqs)
    print("\nNames:", names)
    print("Matrix:")
    for r in mat:
        print(r)
    print("\nASCII visualization:")
    print(ascii_visual(mat))
    if input("\nSave CSV? (y/n): ").lower().startswith("y"):
        save_matrix_csv(names, mat, "hamming_matrix.csv")
        print("Saved hamming_matrix.csv")
